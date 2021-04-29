!$Id:$
!----------------------------------------------------------------------
!
! Code to read city lights file and make a mask
! file in ENVI format, downloaded from:
! http://ngdc.noaa.gov/eog/viirs/download_viirs_ntl.html
!
! Denis Botambekov  -  denis.botambekov@ssec.wisc.edu
!
! April, 2020
!
!----------------------------------------------------------------------

module CITY_MASK_MOD


!--- use statements
use PIXEL_COMMON_MOD, only: & 
      Ch &
    , Sfc &
    , Nav &
    , Geo &
    , Ancil_Data_Dir

use CONSTANTS_MOD, only: &
      int8 &
    , int1 &
    , sym &
    , EXE_PROMPT

use CX_NETCDF4_MOD,only: &
      OPEN_NETCDF &
    , CLOSE_NETCDF &
    , READ_NETCDF 

use NUMERICAL_ROUTINES_MOD,only: &
      FIND_BOUNDS

implicit none

private
public :: READ_CITY_MASK


contains

!----------------------------------------------------------------------
subroutine READ_CITY_MASK()

character(len=1020) :: City_Lights_File
integer :: Ncid
integer :: Nx, Ny, i, j, Temp
integer :: Ilat1, Ilat2, Ilon1, Ilon2, Ilat, Ilon, Ilat_Ad, Ilon_Ad, &
             Ilon1_2, Ilon2_2, Num_Lat, Num_Lon
integer(kind=int1) :: Dateline_Flg
integer, dimension(2) :: Start_2d, Stride_2d, Edge_2d, &
                         Start_2d_2, Stride_2d_2, Edge_2d_2
real, dimension(:,:), allocatable :: CityLights_Read, CityLights_Read2
real :: First_Lat, First_Lon, Res
real :: Wlon, Elon, Slat, Nlat
real, parameter :: CITY_LIGHTS_THRESH = 2.5e-08

! --- original file
!/apollo/cloud/Ancil_Data/clavrx_ancil_data/static/dnb_ancils/
!20150101-20151231_full_globe_vcmcfg_v10.avg_rade9.dat

! --- initialize
Sfc%City_Mask = 0
Num_Lat = 33600
Num_Lon = 86400
First_Lat =   75.
First_Lon = -180.
Res  = 15 * 1./3600
Nx = size(Nav%Lat,1)
Ny = size(Nav%Lat,2)


! --- location of the input file
City_Lights_File = trim(Ancil_Data_Dir)//'/static/dnb_ancils/city_lights.nc'

! --- find bounds
call FIND_BOUNDS(Nav%Lat,Nav%Lon,Wlon,Elon,Slat,Nlat,Dateline_Flg)

! --- if granule north of First_Lat exit
if (Slat > First_Lat) return

!--- open file
call OPEN_NETCDF(trim(City_Lights_File), Ncid)

! --- if granule not crossing dateline
if (dateline_flg == 0) then

   ! --- find slice location
   Ilat1 = max(0,min(Num_Lat,int(abs(Nlat - First_Lat)/Res) + 0))
   Ilat2 = max(0,min(Num_Lat,int(abs(Slat - First_Lat)/Res) + 0))

   Ilon1 = max(0,min(Num_Lon,int(abs(Wlon - First_Lon)/Res) + 0))
   Ilon2 = max(0,min(Num_Lon,int(abs(Elon - First_Lon)/Res) + 0))

   if (Ilat1 > Ilat2) then
      Temp = Ilat1
      Ilat1 = Ilat2
      Ilat2 = Temp
   end if

   if (Ilon1 > Ilon2) then
      Temp = Ilon1
      Ilon1 = Ilon2
      Ilon2 = Temp
   end if

   Start_2d = (/Ilon1, Ilat1/)
   Stride_2d = (/1, 1/)
   Edge_2d = (/(Ilon2-Ilon1)+1, (Ilat2-Ilat1)+1/)

   ! --- allocate buffer
   allocate(CityLights_Read(Edge_2d(1), Edge_2d(2)))

   ! --- call read segment
   call READ_NETCDF(Ncid, Start_2d, Stride_2d, Edge_2d, 'city_lights', CityLights_Read)


   ! --- scale to sensor resolution
   do j = 1, Ny
     do i = 1, Nx

        if (Geo%Space_Mask(i,j) == sym%NO_SPACE) then

           Ilat = max(1,min(Num_Lat,int(abs(Nav%Lat(i,j) - First_Lat)/Res) + 1))
           Ilon = max(1,min(Num_Lon,int(abs(Nav%Lon(i,j) - First_Lon)/Res) + 1))
           Ilat_ad = max(1,min((Ilat - Start_2d(2)) + 1,size(CityLights_Read,2)))
           Ilon_ad = max(1,min((Ilon - Start_2d(1)) + 1,size(CityLights_Read,1)))
           if (CityLights_Read(Ilon_ad,Ilat_ad) > CITY_LIGHTS_THRESH) &
                 Sfc%City_Mask(i,j) = 1

         end if

      end do
   end do

   ! --- deallocate
   if (allocated(CityLights_Read)) deallocate(CityLights_Read)

else ! dateline flag ne 0

   ! --- find slice location
   Ilat1 = max(0,min(Num_Lat,int(abs(Nlat - First_Lat)/Res) + 0))
   Ilat2 = max(0,min(Num_Lat,int(abs(Slat - First_Lat)/Res) + 0))

   Ilon1 = max(0,min(Num_Lon,int(abs(Wlon - First_Lon)/Res) + 0))
   Ilon2 = max(0,min(Num_Lon,int(abs(180.0 - First_Lon)/Res) + 0))

   Ilon1_2 = max(0,min(Num_Lon,int(abs(-180.0 - First_Lon)/Res) + 0))
   Ilon2_2 = max(0,min(Num_Lon,int(abs((Elon-360.0) - First_Lon)/Res) + 0))

   if (Ilat1 > Ilat2) then
      Temp = Ilat1
      Ilat1 = Ilat2
      Ilat2 = Temp
   end if

   if (Ilon1 > Ilon2) then
      Temp = Ilon1
      Ilon1 = Ilon2
      Ilon2 = Temp
   end if

   if (Ilon1_2 > Ilon2_2) then
      Temp = Ilon1_2
      Ilon1_2 = Ilon2_2
      Ilon2_2 = Temp
   end if

   Start_2d = (/Ilon1, Ilat1/)
   Stride_2d = (/1, 1/)
   Edge_2d = (/(Ilon2-Ilon1)+1, (Ilat2-Ilat1)+1/)

   if (Start_2d(1) == 0) Start_2d(1) = 1
   if (Start_2d(2) == 0) Start_2d(2) = 1

   ! --- allocate 1st buffer
   allocate(CityLights_Read(Edge_2d(1), Edge_2d(2)))

   ! --- call read 1st segment
   call READ_NETCDF(Ncid, Start_2d, Stride_2d, Edge_2d, 'city_lights', CityLights_Read)

   Start_2d_2 = (/Ilon1_2, Ilat1/)
   Stride_2d_2 = (/1, 1/)
   Edge_2d_2 = (/(Ilon2_2-Ilon1_2)+1, (Ilat2-Ilat1)+1/)

   if (Start_2d_2(1) == 0) Start_2d_2(1) = 1
   if (Start_2d_2(2) == 0) Start_2d_2(2) = 1

   ! --- allocate 2nd buffer
   allocate(CityLights_Read2(Edge_2d_2(1), Edge_2d_2(2)))

   ! --- call read 2nd segment
   call READ_NETCDF(Ncid, Start_2d_2, Stride_2d_2, Edge_2d_2, 'city_lights', CityLights_Read2)

   ! --- scale to sensor resolution
   do j = 1, Ny
     do i = 1, Nx

        if (Geo%Space_Mask(i,j) == sym%NO_SPACE) then

           Ilat = max(1,min(Num_Lat,int(abs(Nav%Lat(i,j) - First_Lat)/Res) + 1))
           Ilon = max(1,min(Num_Lon,int(abs(Nav%Lon(i,j) - First_Lon)/Res) + 1))

           if (Nav%Lon(i,j) >= 0.0) then
              Ilat_ad = max(1,min((Ilat - Start_2d(2)) + 1,size(CityLights_Read,2)))
              Ilon_ad = max(1,min((Ilon - Start_2d(1)) + 1,size(CityLights_Read,1)))
              if (CityLights_Read(Ilon_ad,Ilat_ad) > CITY_LIGHTS_THRESH) &
                   Sfc%City_Mask(i,j) = 1
           else
              Ilat_ad = max(1,min((Ilat - Start_2d_2(2)) + 1,size(CityLights_Read2,2)))
              Ilon_ad = max(1,min((Ilon - Start_2d_2(1)) + 1,size(CityLights_Read2,1)))
              if (CityLights_Read2(Ilon_ad,Ilat_ad) > CITY_LIGHTS_THRESH) &
                   Sfc%City_Mask(i,j) = 1
           end if

         end if

      end do
   end do

   ! --- deallocate
   if (allocated(CityLights_Read)) deallocate(CityLights_Read)
   if (allocated(CityLights_Read2)) deallocate(CityLights_Read2)

endif

! --- close file
call CLOSE_NETCDF(Ncid)


end subroutine READ_CITY_MASK
!----------------------------------------------------------------------

end module CITY_MASK_MOD

