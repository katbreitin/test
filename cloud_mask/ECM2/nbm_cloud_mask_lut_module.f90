!$Id: nbm_cloud_mask_lut_module.f90 
!----------------------------------------------------------------------
! MODULE name: NBM_CLOUD_MASK_LUT_MODULE
! 
! Routines for the naive Bayesian cloud mask LUT reading
! Version 2.0
!
! Authors: Andrew Heidinger, NOAA/NESDIS
!          Andi Walther, CIMSS
!          Denis Botambekov, CIMSS
!          William Straka, CIMSS
!
! DEPENDENCIES: Services_Module, NetCDF_read_Module, File_Tools
!
! SIDE EFFECTS: None
!
!-------------------------------------------------------------------

module NBM_CLOUD_MASK_LUT_MODULE


 use NB_CLOUD_MASK_SERVICES

 use NB_CLOUD_MASK_NETCDF_READ_MODULE, only: &
     OPEN_NETCDF &
   , CLOSE_NETCDF &
   , READ_NETCDF_2D_CHAR &
   , READ_NETCDF_1D_INT  &
   , READ_NETCDF_1D_REAL &
   , READ_NETCDF_2D_REAL &
   , READ_NETCDF_3D_REAL &
   , READ_NETCDF_4D_REAL &
   , GET_GROUP_ID &
   , READ_NETCDF_GLOBAL_ATTRIBUTE_I4 &
   , READ_NETCDF_GLOBAL_ATTRIBUTE_R4 &
   , READ_NETCDF_GLOBAL_ATTRIBUTE_CHAR &
   , READ_NETCDF_DIMENSION &
   , DEFAULT_NAME_LENGTH

 use FILE_TOOLS, only: &
     FILE_TEST

 implicit none

 public NBM_CLOUD_MASK_COMPUTE_PRIOR
 public NBM_CLOUD_MASK_LUT_READ
 public RESET_NBM_CLOUD_MASK_LUT
 public RESET_NBM_CLOUD_MASK_PRIOR_LUT

 logical, public, save:: Is_Classifiers_Read = .false.
 logical, public, save:: Is_Prior_Read_M = .false.
! This has to be here to be available to bridge, science code. - WCS3
 type(Classifier), dimension(:), allocatable, public, save :: Lut
 type(Mask_Threshold), public, save :: Mask_Thresh

 contains

!-------------------------------------------------------------------

subroutine NBM_CLOUD_MASK_COMPUTE_PRIOR(Prior_File_Name, Lon, Lat, Month, Prior_Probability)

character(len=*), intent(in) :: Prior_File_Name
real, dimension(:,:), intent(in) :: Lon
real, dimension(:,:), intent(in) :: Lat
integer(kind=int2), intent(in) :: Month
real(kind=real4), dimension(:,:), intent(inout) :: Prior_Probability

real :: Nlon_Prior, Nlat_Prior, Nmonths_Prior, Ndiurnal_Prior
real :: Dlon_Prior, Dlat_Prior
real :: Lon_Min_Prior, Lon_Max_Prior
real :: Lat_Min_Prior, Lat_Max_Prior
real(kind=real4), dimension(:,:,:,:), allocatable :: Prior_Table
integer :: Ncid
integer :: Nx, Ny
integer, parameter:: Idiurnal = 1    !1 = daily averaged
integer :: i,j,ilon,ilat
integer, dimension(4) :: Dims
logical :: File_Exist

! --- initiate
Is_Prior_Read_M = .false.
File_Exist = .false.

!--- check if file exists
File_Exist = FILE_TEST(trim(Prior_File_Name))
if (.not. File_Exist) then
  print *,'PRIOR CLOUD MASK LUT IS NOT FOUND, STOPPING'
  stop
endif

!--- open file
call OPEN_NETCDF(trim(Prior_File_Name), Ncid)

!--- read attributes
call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Ncid, "number_longitudes", Nlon_Prior)
call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Ncid, "number_latitudes", Nlat_Prior)
call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Ncid, "number_months", Nmonths_Prior)
call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Ncid, "number_times", Ndiurnal_Prior)
call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Ncid, "longitude_spacing", Dlon_Prior)
call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Ncid, "latitude_spacing", Dlat_Prior)
call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Ncid, "longitude_min", Lon_Min_Prior)
call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Ncid, "longitude_max", Lon_Max_Prior)
call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Ncid, "latitude_min", Lat_Min_Prior)
call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Ncid, "latitude_max", Lat_Max_Prior)


! --- allocate prior table
Dims = (/int(Nlon_Prior),int(Nlat_Prior),int(Nmonths_Prior),int(Ndiurnal_Prior)/)
allocate(Prior_Table(Dims(1),Dims(2),Dims(3),Dims(4)))

! --- read table
call READ_NETCDF_4D_REAL(Ncid, (/1,1,1,1/), (/1,1,1,1/), Dims, &
               !"cloud_fraction_table", Prior_Table)
               "cloud_fraction_table_smoothed", Prior_Table)

! --- close file
call CLOSE_NETCDF(Ncid)


! --- get dimensions 
Nx = size(Lon,1)
Ny = size(Lon,2)

! --- calculate prior
do i = 1, Nx
   do j = 1, Ny

      ilon = min(int(Nlon_Prior),max(1,int((Lon(i,j) - Lon_Min_Prior) / (Dlon_Prior))+1))
      ilat = min(int(Nlat_Prior),max(1,int((Lat(i,j) - Lat_Min_Prior) / (Dlat_Prior))+1))

      Prior_Probability(i,j) = Prior_Table(Ilon, Ilat, Month, Idiurnal)

   enddo
enddo

Is_Prior_Read_M = .true.

! --- deallocate prior table
if (allocated(Prior_Table)) deallocate(Prior_Table)

end subroutine NBM_CLOUD_MASK_COMPUTE_PRIOR

!-------------------------------------------------------------------

subroutine NBM_CLOUD_MASK_LUT_READ(Lut_File_Full_Path, N_Classifier)

character(len=*), intent(in) :: Lut_File_Full_Path
integer, intent(out) :: N_Classifier

integer :: Ncid, Group_Id
integer :: Class_Idx
integer, dimension(2) :: Start_Read_2d, Count_Read_2d
integer, dimension(3) :: Start_Read_3d, Count_Read_3d
integer, dimension(4) :: Start_Read_4d, Count_Read_4d
real, dimension(:,:), allocatable :: Buffer_2d
real, dimension(:,:,:), allocatable :: Buffer_3d
real, dimension(:,:,:,:), allocatable :: Buffer_4d
real, parameter :: MISSING = -999.0
character(len=100) :: Attr_Name
character(len=100) :: Dim_Name
character(len=100) :: Var_Name
logical :: File_Exist
character(len=100) :: Temp_Name

! --- initiate
Is_Classifiers_Read = .FALSE.
File_Exist = .FALSE.

!--- check if file exists
File_Exist = FILE_TEST(trim(Lut_File_Full_Path))
if (.not. File_Exist) then
  print *,'CLOUD MASK LUT IS NOT FOUND, STOPPING'
  stop
endif

!--- open file
call OPEN_NETCDF(trim(Lut_File_Full_Path), Ncid)

!--- read name length
Dim_Name = 'nlength_class'
call READ_NETCDF_DIMENSION(Ncid, trim(Dim_Name), DEFAULT_NAME_LENGTH)

!--- read number and name of classifiers from file
Attr_Name = 'nclassifiers'
call READ_NETCDF_GLOBAL_ATTRIBUTE_I4(Ncid, trim(Attr_Name), N_Classifier)

! --- read classifier names
Var_Name = 'classifier_names'
allocate(Classifier_Names(N_Classifier),MOLD=Temp_Name(1:DEFAULT_NAME_LENGTH))
call READ_NETCDF_2D_CHAR(Ncid, (/1,1/), Var_Name, Classifier_Names)
!do Class_Idx = 1, N_Classifier
! print *, "Classifier Names = ", Classifier_Names(Class_Idx)
!enddo

!---- allocate main structure based on number of classifiers
allocate(Lut(N_Classifier))

! --- read tut/rut thresh
allocate(Mask_Thresh%Rut_Clear_Prob_Clear_Thresh(7))
allocate(Mask_Thresh%Tut_Clear_Prob_Clear_Thresh(7)) 
call READ_NETCDF_1D_REAL(Ncid, (/1/), (/1/), (/7/), 'rut_clear_prob_clear_thresh', Mask_Thresh%Rut_Clear_Prob_Clear_Thresh)
call READ_NETCDF_1D_REAL(Ncid, (/1/), (/1/), (/7/), 'tut_clear_prob_clear_thresh', Mask_Thresh%Tut_Clear_Prob_Clear_Thresh)


! --- loop over classifiers
do Class_Idx = 1, N_Classifier

   ! --- open group
   call GET_GROUP_ID(Ncid, trim(Classifier_Names(Class_Idx)), Group_Id)

   !allocate names
   allocate(character(len=DEFAULT_NAME_LENGTH) :: Lut(Class_Idx)%Class_Xname)
   allocate(character(len=DEFAULT_NAME_LENGTH) :: Lut(Class_Idx)%Class_Yname)
   allocate(character(len=DEFAULT_NAME_LENGTH) :: Lut(Class_Idx)%Class_Zname)

   ! --- initialize
   Lut(Class_Idx)%Zen_Min = MISSING
   Lut(Class_Idx)%Zen_Max = MISSING
   Lut(Class_Idx)%Solzen_Min = MISSING
   Lut(Class_Idx)%Solzen_Max = MISSING
   Lut(Class_Idx)%Solglintzen_Min = MISSING
   Lut(Class_Idx)%Solglintzen_Max = MISSING
   Lut(Class_Idx)%Zsfc_Min = MISSING
   Lut(Class_Idx)%Zsfc_Max = MISSING
   Lut(Class_Idx)%Zsfc_Std_Min = MISSING
   Lut(Class_Idx)%Zsfc_Std_Max = MISSING
   Lut(Class_Idx)%Tsfc_Min = MISSING
   Lut(Class_Idx)%Tsfc_Max = MISSING
   Lut(Class_Idx)%Tpw_Min = MISSING
   Lut(Class_Idx)%Tpw_Max = MISSING
   Lut(Class_Idx)%Rut_Solzen_Thresh = MISSING

   ! --- read in group attributes
   call READ_NETCDF_GLOBAL_ATTRIBUTE_I4(Group_Id, 'rank', Lut(Class_Idx)%Rank)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_I4(Group_Id, 'nchan_used', Lut(Class_Idx)%Nchan_Used)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_I4(Group_Id, 'nsfc', Lut(Class_Idx)%N_Sfc)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_CHAR(Group_Id, 'x_name', Lut(Class_Idx)%Class_Xname)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_I4(Group_Id, 'nbins_x', Lut(Class_Idx)%Nbins_X)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'x_min', Lut(Class_Idx)%X_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'x_bin', Lut(Class_Idx)%X_Bin)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_CHAR(Group_Id, 'y_name', Lut(Class_Idx)%Class_Yname)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_I4(Group_Id, 'nbins_y', Lut(Class_Idx)%Nbins_Y)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'y_min', Lut(Class_Idx)%Y_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'y_bin', Lut(Class_Idx)%Y_Bin)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_CHAR(Group_Id, 'z_name', Lut(Class_Idx)%Class_Zname)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_I4(Group_Id, 'nbins_z', Lut(Class_Idx)%Nbins_Z)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'z_min', Lut(Class_Idx)%Z_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'z_bin', Lut(Class_Idx)%Z_Bin)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'zenith_minimum', Lut(Class_Idx)%Zen_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'zenith_maximum', Lut(Class_Idx)%Zen_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'solar_zenith_minimum', Lut(Class_Idx)%Solzen_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'solar_zenith_maximum', Lut(Class_Idx)%Solzen_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'solar_glint_zenith_minimum', Lut(Class_Idx)%Solglintzen_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'solar_glint_zenith_maximum', Lut(Class_Idx)%Solglintzen_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'solglint_mask_minimum', Lut(Class_Idx)%Solglint_Mask_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'solglint_mask_maximum', Lut(Class_Idx)%Solglint_Mask_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'surface_elevation_minimum', Lut(Class_Idx)%Zsfc_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'surface_elevation_maximum', Lut(Class_Idx)%Zsfc_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'surface_elevation_stddev_minimum', Lut(Class_Idx)%Zsfc_Std_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'surface_elevation_stddev_maximum', Lut(Class_Idx)%Zsfc_Std_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'surface_temperature_minimum', Lut(Class_Idx)%Tsfc_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'surface_temperature_maximum', Lut(Class_Idx)%Tsfc_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'total_precipitable_water_minimum', Lut(Class_Idx)%Tpw_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'total_precipitable_water_maximum', Lut(Class_Idx)%Tpw_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'snow_class_minimum', Lut(Class_Idx)%Snow_Class_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'snow_class_maximum', Lut(Class_Idx)%Snow_Class_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'lunar_zenith_minimum', Lut(Class_Idx)%Lunzen_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'lunar_zenith_maximum', Lut(Class_Idx)%Lunzen_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'lunar_glint_zenith_minimum', Lut(Class_Idx)%Lunglintzen_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'lunar_glint_zenith_maximum', Lut(Class_Idx)%Lunglintzen_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'lunglint_mask_minimum', Lut(Class_Idx)%Lunglint_Mask_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'lunglint_mask_maximum', Lut(Class_Idx)%Lunglint_Mask_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'solar_scattering_angle_minimum', Lut(Class_Idx)%Solscatang_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'solar_scattering_angle_maximum', Lut(Class_Idx)%Solscatang_Max)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'coast_mask_minimum', Lut(Class_Idx)%Coast_Mask_Min)
   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'coast_mask_maximum', Lut(Class_Idx)%Coast_Mask_Max)

   ! --- temporary set before all LUTs will have these!!! Denis B. 2021-01-13
   !call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'city_mask_minimum', Lut(Class_Idx)%City_Mask_Min)
   !call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'city_mask_maximum', Lut(Class_Idx)%City_Mask_Max)
   Lut(Class_Idx)%City_Mask_Min = 0
   Lut(Class_Idx)%City_Mask_Max = 1
   !call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'moon_illum_frac_minimum', Lut(Class_Idx)%Moon_Illum_Frac_Min)
   !call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'moon_illum_frac_maximum', Lut(Class_Idx)%Moon_Illum_Frac_Max)
   Lut(Class_Idx)%Moon_Illum_Frac_Min = 0.
   Lut(Class_Idx)%Moon_Illum_Frac_Max = 1.

   call READ_NETCDF_GLOBAL_ATTRIBUTE_R4(Group_Id, 'rut_solzen_thresh', Lut(Class_Idx)%Rut_Solzen_Thresh) 

   ! --- read wave length
   allocate(Lut(Class_Idx)%Wvl(Lut(Class_Idx)%Nchan_Used))
   call READ_NETCDF_1D_INT(Group_Id, (/1/), (/1/), (/Lut(Class_Idx)%Nchan_Used/), 'channel_wvl', Lut(Class_Idx)%Wvl)

   ! --- read on_flag
   allocate(Lut(Class_Idx)%On_Flag(Lut(Class_Idx)%N_Sfc))
   call READ_NETCDF_1D_INT(Group_Id, (/1/), (/1/), (/Lut(Class_Idx)%N_Sfc/), 'on_flag', Lut(Class_Idx)%On_Flag)
   
   !alloate and read in prior LUT - WCS3
   allocate(Lut(Class_Idx)%Cloud_fraction(Lut(Class_Idx)%N_Sfc))
   call READ_NETCDF_1D_REAL(Group_Id, (/1/), (/1/), (/Lut(Class_Idx)%N_Sfc/), 'cloud_fraction', &
                             Lut(Class_Idx)%Cloud_fraction)

   allocate(Lut(Class_Idx)%Ice_Fraction(Lut(Class_Idx)%N_Sfc))
   call READ_NETCDF_1D_REAL(Group_Id, (/1/), (/1/), (/Lut(Class_Idx)%N_Sfc/), 'ice_fraction', &
                             Lut(Class_Idx)%Ice_Fraction)

   allocate(Lut(Class_Idx)%Water_Fraction(Lut(Class_Idx)%N_Sfc))
   call READ_NETCDF_1D_REAL(Group_Id, (/1/), (/1/), (/Lut(Class_Idx)%N_Sfc/), 'water_fraction', &
                             Lut(Class_Idx)%Water_Fraction)

   
   ! --- allocatable based on sizes in attributes
   if (Lut(Class_Idx)%Rank == 1) then
      allocate(Lut(Class_Idx)%Clear_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%N_Sfc,1,1))
      allocate(Lut(Class_Idx)%Water_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%N_Sfc,1,1))
      allocate(Lut(Class_Idx)%Ice_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%N_Sfc,1,1))
      allocate(Lut(Class_Idx)%Obs_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%N_Sfc,1,1))
      allocate(Buffer_2d(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%N_Sfc))
      ! --- read in lut tables
      Start_Read_2d = 1
      Count_Read_2d(1) = Lut(Class_Idx)%Nbins_X
      Count_Read_2d(2) = Lut(Class_Idx)%N_Sfc
      call READ_NETCDF_2D_REAL(Group_Id, Start_Read_2d, (/1,1/), Count_Read_2d, 'clear_table', Buffer_2d)
      Lut(Class_Idx)%Clear_Table(:,:,1,1) = Buffer_2d(:,:)

      call READ_NETCDF_2D_REAL(Group_Id, Start_Read_2d, (/1,1/), Count_Read_2d, 'water_table', Buffer_2d)
      Lut(Class_Idx)%Water_Table(:,:,1,1) = Buffer_2d(:,:)

      call READ_NETCDF_2D_REAL(Group_Id, Start_Read_2d, (/1,1/), Count_Read_2d, 'ice_table', Buffer_2d)
      Lut(Class_Idx)%Ice_Table(:,:,1,1) = Buffer_2d(:,:)

      call READ_NETCDF_2D_REAL(Group_Id, Start_Read_2d, (/1,1/), Count_Read_2d, 'obs_table', Buffer_2d)
      Lut(Class_Idx)%Obs_Table(:,:,1,1) = Buffer_2d(:,:)


   elseif (Lut(Class_Idx)%Rank == 2) then
      allocate(Lut(Class_Idx)%Clear_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%Nbins_Y, &
                                        Lut(Class_Idx)%N_Sfc,1))
      allocate(Lut(Class_Idx)%Water_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%Nbins_Y, &
                                        Lut(Class_Idx)%N_Sfc,1))
      allocate(Lut(Class_Idx)%Ice_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%Nbins_Y, &
                                        Lut(Class_Idx)%N_Sfc,1))
      allocate(Lut(Class_Idx)%Obs_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%Nbins_Y, &
                                        Lut(Class_Idx)%N_Sfc,1))
      allocate(Buffer_3d(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%Nbins_Y,Lut(Class_Idx)%N_Sfc))
      ! --- read in lut tables
      Start_Read_3d = 1
      Count_Read_3d(1) = Lut(Class_Idx)%Nbins_X
      Count_Read_3d(2) = Lut(Class_Idx)%Nbins_Y
      Count_Read_3d(3) = Lut(Class_Idx)%N_Sfc
      call READ_NETCDF_3D_REAL(Group_Id, Start_Read_3d, (/1,1,1/), Count_Read_3d, 'clear_table', Buffer_3d)
      Lut(Class_Idx)%Clear_Table(:,:,:,1) = Buffer_3d(:,:,:)

      call READ_NETCDF_3D_REAL(Group_Id, Start_Read_3d, (/1,1,1/), Count_Read_3d, 'water_table', Buffer_3d)
      Lut(Class_Idx)%Water_Table(:,:,:,1) = Buffer_3d(:,:,:)

      call READ_NETCDF_3D_REAL(Group_Id, Start_Read_3d, (/1,1,1/), Count_Read_3d, 'ice_table', Buffer_3d)
      Lut(Class_Idx)%Ice_Table(:,:,:,1) = Buffer_3d(:,:,:)

      call READ_NETCDF_3D_REAL(Group_Id, Start_Read_3d, (/1,1,1/), Count_Read_3d, 'obs_table', Buffer_3d)
      Lut(Class_Idx)%Obs_Table(:,:,:,1) = Buffer_3d(:,:,:)


   elseif (Lut(Class_Idx)%Rank == 3) then
      allocate(Lut(Class_Idx)%Clear_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%Nbins_Y, &
                                        Lut(Class_Idx)%Nbins_Z,Lut(Class_Idx)%N_Sfc))
      allocate(Lut(Class_Idx)%Water_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%Nbins_Y, &
                                        Lut(Class_Idx)%Nbins_Z,Lut(Class_Idx)%N_Sfc))
      allocate(Lut(Class_Idx)%Ice_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%Nbins_Y, &
                                        Lut(Class_Idx)%Nbins_Z,Lut(Class_Idx)%N_Sfc))
      allocate(Lut(Class_Idx)%Obs_Table(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%Nbins_Y, &
                                        Lut(Class_Idx)%Nbins_Z,Lut(Class_Idx)%N_Sfc))
      allocate(Buffer_4d(Lut(Class_Idx)%Nbins_X,Lut(Class_Idx)%Nbins_Y, &
                                        Lut(Class_Idx)%Nbins_Z, Lut(Class_Idx)%N_Sfc))
      ! --- read in lut tables
      Start_Read_4d = 1
      Count_Read_4d(1) = Lut(Class_Idx)%Nbins_X
      Count_Read_4d(2) = Lut(Class_Idx)%Nbins_Y
      Count_Read_4d(3) = Lut(Class_Idx)%Nbins_Z
      Count_Read_4d(4) = Lut(Class_Idx)%N_Sfc
      call READ_NETCDF_4D_REAL(Group_Id, Start_Read_4d, (/1,1,1,1/), Count_Read_4d, 'clear_table', Buffer_4d)
      Lut(Class_Idx)%Clear_Table(:,:,:,:) = Buffer_4d(:,:,:,:)
      call READ_NETCDF_4D_REAL(Group_Id, Start_Read_4d, (/1,1,1,1/), Count_Read_4d, 'water_table', Buffer_4d)
      Lut(Class_Idx)%Water_Table(:,:,:,:) = Buffer_4d(:,:,:,:)
      call READ_NETCDF_4D_REAL(Group_Id, Start_Read_4d, (/1,1,1,1/), Count_Read_4d, 'ice_table', Buffer_4d)
      Lut(Class_Idx)%Ice_Table(:,:,:,:) = Buffer_4d(:,:,:,:)
      call READ_NETCDF_4D_REAL(Group_Id, Start_Read_4d, (/1,1,1,1/), Count_Read_4d, 'obs_table', Buffer_4d)
      Lut(Class_Idx)%Obs_Table(:,:,:,:) = Buffer_4d(:,:,:,:)

   else
      print *,'Rank is > 3, Stopping'
      stop
   endif

   ! --- deallocate buffer
   if (allocated(Buffer_2d)) deallocate(Buffer_2d)
   if (allocated(Buffer_3d)) deallocate(Buffer_3d)
   if (allocated(Buffer_4d)) deallocate(Buffer_4d)


enddo ! loop over n class


! --- allocate prob/conf clear/cloudy thresholes 
allocate(Mask_Thresh%Conf_Clear_Prob_Clear_Thresh(Lut(1)%N_Sfc))
allocate(Mask_Thresh%Prob_Clear_Prob_Cloudy_Thresh(Lut(1)%N_Sfc))
allocate(Mask_Thresh%Prob_Cloudy_Conf_Cloudy_Thresh(Lut(1)%N_Sfc))

! --- set to missing
Mask_Thresh%Conf_Clear_Prob_Clear_Thresh = MISSING
Mask_Thresh%Prob_Clear_Prob_Cloudy_Thresh = MISSING
Mask_Thresh%Prob_Cloudy_Conf_Cloudy_Thresh = MISSING

! --- read prob/conf clear/cloudy thresholes 
call READ_NETCDF_1D_REAL(Ncid, (/1/), (/1/),(/Lut(1)%N_Sfc/), &
                'conf_clear_prob_clear_thresh', Mask_Thresh%Conf_Clear_Prob_Clear_Thresh)
call READ_NETCDF_1D_REAL(Ncid, (/1/), (/1/),(/Lut(1)%N_Sfc/), &
                'prob_clear_prob_cloud_thresh', Mask_Thresh%Prob_Clear_Prob_Cloudy_Thresh)
call READ_NETCDF_1D_REAL(Ncid, (/1/), (/1/),(/Lut(1)%N_Sfc/), &
                'prob_cloud_conf_cloud_thresh', Mask_Thresh%Prob_Cloudy_Conf_Cloudy_Thresh)

!

!--- close file
call CLOSE_NETCDF(Ncid)

Is_Classifiers_Read = .true.

!deallocate(Classifier_Names)
if (allocated(Buffer_2d)) deallocate(Buffer_2d)
if (allocated(Buffer_2d)) deallocate(Buffer_3d)
if (allocated(Buffer_2d)) deallocate(Buffer_4d)

end subroutine NBM_CLOUD_MASK_LUT_READ

!-----------------------------------------------------------------------------------------
! This routine deallocates all LUT arrays and resets Is_Classifiers_Read to be false
!-----------------------------------------------------------------------------------------
subroutine RESET_NBM_CLOUD_MASK_LUT()

deallocate(Classifier_Names)
deallocate(Mask_Thresh%Conf_Clear_Prob_Clear_Thresh)
deallocate(Mask_Thresh%Prob_Clear_Prob_Cloudy_Thresh)
deallocate(Mask_Thresh%Prob_Cloudy_Conf_Cloudy_Thresh)
deallocate(Mask_Thresh%Rut_Clear_Prob_Clear_Thresh)
deallocate(Mask_Thresh%Tut_Clear_Prob_Clear_Thresh)
deallocate(Lut)

Is_Classifiers_Read = .false.

end subroutine RESET_NBM_CLOUD_MASK_LUT

!-----------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------
subroutine RESET_NBM_CLOUD_MASK_PRIOR_LUT()

 Is_Prior_Read_M = .false.

end subroutine RESET_NBM_CLOUD_MASK_PRIOR_LUT


!=========================================================================================

end module NBM_CLOUD_MASK_LUT_MODULE

