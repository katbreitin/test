! $Id: acha_shadow_module.f90 3083 2018-12-17 18:36:41Z mfoster $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: acha_shadow_module.f90 (src)
!       ACHA_SHADOW (program)
!
! PURPOSE: this module computes the geometrical shadow routine
!
! DESCRIPTION: 
!
! AUTHORS:  Andi Walther, UW/SSEC/CIMSS
!
! COPYRIGHT
! THIS SOFTWARE AND ITS DOCUMENTATION ARE CONSIDERED TO BE IN THE PUBLIC
! DOMAIN AND THUS ARE AVAILABLE FOR UNRESTRICTED PUBLIC USE. THEY ARE
! FURNISHED "AS IS." THE AUTHORS, THE UNITED STATES GOVERNMENT, ITS
! INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND AGENTS MAKE NO WARRANTY,
! EXPRESS OR IMPLIED, AS TO THE USEFULNESS OF THE SOFTWARE AND
! DOCUMENTATION FOR ANY PURPOSE. THEY ASSUME NO RESPONSIBILITY (1) FOR
! THE USE OF THE SOFTWARE AND DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL
! SUPPORT TO USERS.
!
! Public routines used in this MODULE:
!
!--------------------------------------------------------------------------------------
module ACHA_SHADOW

 use ACHA_NUM_MOD, only: FINDGEN

 implicit none

 public:: CLOUD_SHADOW_RETR

 contains

!==============================================================================
! Cloud Shadow Routine
!
! -1 = shadow estimation not possible
!  0 = no shadow
!  1 = shadow determined from geometrical analysis
!  2 = shadow determined from spectral analysis
!  3 = shadow determined from both geometrical and spectral analysis
!
!==============================================================================
subroutine CLOUD_SHADOW_RETR (  &
           Zt, &
           Zb, &
           Zs, &
           SolAz, &
           SolZen, &
           Lat, &
           Lon, &
           CldMask, &
           CldOpd, &
           Cloud_Shadow )
      
      implicit none 
      real, intent(in), dimension(:,:) :: Zt
      real, intent(in), dimension(:,:) :: Zb
      real, intent(in), dimension(:,:) :: Zs
      real, intent(in), dimension(:,:) :: Solaz
      real, intent(in), dimension(:,:) :: SolZen
      real, intent(in), dimension(:,:) :: Lat
      real, intent(in), dimension(:,:) :: Lon
      integer(1), intent(in), dimension(:,:) :: CldMask
      real, intent(in), dimension(:,:) :: CldOpd
      integer(1), intent(out), dimension(:,:) :: Cloud_Shadow
      
      real, parameter :: PI = 3.1415926535897
      real, parameter :: DTOR = PI/180.
      real :: Lon_Spacing_Per_m
      real, parameter:: LAT_SPACING_PER_M = 8.9932e-06   ! ( = 1.0/111000.0 m )
      real, parameter:: MISSING_VALUE_REAL4 = -999.0
      integer, parameter:: Npts_Max = 50
      
      integer :: i,j
      real :: Delta_Lon, Delta_Lat, Delta_Lon_Per_Pixel, Delta_Lat_Per_Pixel
      real :: Total_Displacement, Deg_Dist_per_m
      integer :: i_dim, j_dim
      real :: x_start, x_cloud, x_end
      real :: y_start, y_cloud, y_end
      real :: temp, z_start_end, Pixel_Width
      integer:: Nx, Ny, Npts
      real :: Lat_Start, Lat_End, Lon_Start, Lon_End
      real, dimension(Npts_Max):: Z_Temp, X_Temp, Y_Temp
      integer, dimension(Npts_Max):: I_Temp, J_Temp
      
      !--- initialize output
      Cloud_Shadow = 0

      Nx = size(lat,1)
      Ny = size(lat,2)

      do i = 2,Nx-1

        x_cloud = float(i)

        do j = 2,Ny-1
      
          y_cloud = float(j)

          Delta_Lon_Per_Pixel = Lon(i+1,j) - Lon(i,j)
          Delta_Lat_Per_Pixel = Lat(i,j+1) - Lat(i,j)
          Pixel_Width = 0.5*abs(Delta_Lat_Per_Pixel + Delta_Lon_Per_Pixel)

          !--- exclude clear
          if (CldMask(i,j) == 0) cycle
          if (CldMask(i,j) == 1) cycle

          !-- exclude missing
          if (Zt(i,j) == MISSING_VALUE_REAL4) cycle

          Lon_Spacing_Per_m = Lat_Spacing_Per_m / cos(Lat(i,j)*DTOR)
          Deg_Dist_per_m = 0.5*(Lat_Spacing_Per_m + Lon_Spacing_Per_m)

          !--- find end point of shadow from cloud top
          temp = tan(Solzen(i,j)*DTOR)* (Zt(i,j) - Zs(i,j))
          Total_Displacement = max(0.0,temp)
          Delta_Lon = sin(Solaz(i,j)*DTOR)*Total_Displacement * Lon_Spacing_Per_m
          Delta_Lat = cos(Solaz(i,j)*DTOR)*Total_Displacement * Lat_Spacing_Per_m
          Lon_End = Lon(i,j) + Delta_Lon
          Lat_End = Lat(i,j) + Delta_Lat

          x_end = float(min(Nx,max(1,int(i + Delta_Lon / Delta_Lon_Per_Pixel))))
          y_end = float(min(Ny,max(1,int(j + Delta_Lat / Delta_Lat_Per_Pixel))))

          !--- find start point of shadow from cloud base
          Delta_Lon = 0.0
          Delta_Lat = 0.0
          if (Zb(i,j) > Zs(i,j) .and. Zb(i,j) < Zt(i,j)) then
             Total_Displacement = max(0.0,tan(Solzen(i,j)*DTOR)* (Zb(i,j) - Zs(i,j)))
             Delta_Lon = sin(Solaz(i,j)*DTOR)*Total_Displacement * Lon_Spacing_Per_m
             Delta_Lat = cos(Solaz(i,j)*DTOR)*Total_Displacement * Lat_Spacing_Per_m
          endif
          Lon_Start = Lon(i,j) + Delta_Lon
          Lat_Start = Lat(i,j) + Delta_Lat
          x_start = float(min(Nx,max(1,int(i + Delta_Lon / Delta_Lon_Per_Pixel))))
          y_start = float(min(Ny,max(1,int(j + Delta_Lat / Delta_Lat_Per_Pixel))))

          if (abs(x_end-x_start) < 1 .and. abs(y_end - y_start) < 1) then
             if (abs(x_cloud - x_start) < 1 .and. abs(y_cloud - y_start) < 1) then
               !print, 'no shadow'
               cycle
             endif
          endif

          !--- find pixels on this line
          Z_start_end = sqrt((X_end - X_start)**2 + (Y_end - Y_start)**2)
          Npts = min(Npts_Max,max(2,2*int(Z_Start_End)))
          Z_temp(1:Npts) = FINDGEN(Npts)*Z_Start_end/float(Npts-1)

          !if (z_start_end gt 50) then stop

          X_Temp(1:Npts) = X_start + sin(Solaz(i,j)*DTOR) * Z_temp(1:Npts)
          where(X_Temp < 1)
                X_Temp = 1
          endwhere
          where(X_Temp > Nx)
                X_Temp = Nx
          endwhere

          Y_Temp(1:Npts) = X_start + cos(Solaz(i,j)*DTOR) * Z_temp(1:Npts)
          where(Y_Temp < 1)
                Y_Temp = 1
          endwhere
          where(Y_Temp > Ny)
                Y_Temp = Ny
          endwhere

          i_temp(1:Npts) = nint(X_Temp(1:Npts))
          j_temp(1:Npts) = nint(Y_Temp(1:Npts))

          Cloud_Shadow(i_temp(1:Npts),j_temp(1:Npts)) = 1


        enddo

      enddo


   end subroutine CLOUD_SHADOW_RETR

!-----------------------------------------------------------
! end of module
!-----------------------------------------------------------
end module ACHA_SHADOW
