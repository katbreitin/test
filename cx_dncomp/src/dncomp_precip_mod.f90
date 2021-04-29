! $Id:$
! CREATED Andi Walther {date} 
module dncomp_precip_mod
   implicit none
contains

   subroutine compute_precip ( cod, cps, rain_column, rain_prob, rain_rate )
      real, intent(in) :: cod
      real, intent(in) :: cps
      real, intent(in) :: rain_column
      real, intent(out) :: rain_prob
      real, intent(out) :: rain_rate

      real, parameter :: dH_0 = 0.6 ! km
      real, parameter :: CWP_0 = 120.0 ! g/m^2
      real, parameter:: Lapse_Rate = 6.5 !K/km
      real, parameter:: Alpha = 1.6 !dimensionless
      real, parameter:: C = 1.0 !mm / hour
      real, parameter:: CWP_T = 150.0 !g/m^2
      real, parameter:: Ceps_T = 15.0 !micron
      real, parameter:: dh_max = 7.0
      real, parameter:: COD_T = 1.

      real :: rain_rate_max


      real :: cwp
      rain_prob = 0.
      rain_rate = 0.
      cwp = cod * cps * 5./9.
      if ( cwp .lt. CWP_T .or.  CPS .lt. CEPS_T .or. COD .lt. COD_T ) return


      rain_rate = (((cwp - cwp_0)/cwp_0 ) ** alpha ) / rain_column
      Rain_Rate_Max = 5.0 + rain_column**Alpha
      rain_rate = min ( rain_rate, rain_rate_max)

      rain_prob = 1.

   end subroutine

   subroutine compute_rain_column ( ctt, rain_column)
      real,intent(in) :: ctt(:,:)
      real,intent(out) :: rain_column(:,:) ! km
      real, parameter:: Lapse_Rate = 6.5 !K/km
      integer :: nx
      integer :: ny
      integer :: size_arr(2)
      integer :: y_idx
      integer :: x_idx

      integer:: x_Idx_1
      integer:: x_Idx_2

      integer:: y_Idx_1
      integer:: y_Idx_2

      real, parameter:: dh_max = 7.0
      real :: ctt_max
      real, parameter :: dH_0 = 0.6 ! km
      integer, parameter:: N_box = 50 ! may have to be defined according resolution

      size_arr = shape ( ctt)
      nx = size_arr(1)
      ny = size_arr(2)

      ! compute rain_column
      ! first find maximal ctt in area
      y_loop: do y_idx = 1, ny
         x_loop: do x_idx = 1, nx
            if ( ctt (  x_idx,y_idx ) .lt. 0 ) cycle
                  !--- define box to look for maximum cloud temperature
            x_Idx_1  = min(nx-1,max(1,x_Idx - N_box /2))
            x_Idx_2  = min(nx,max(2,x_Idx + N_box /2))
            y_Idx_1  = min(ny-1,max(1,y_Idx - N_box /2))
            y_Idx_2  = min(ny,max(2,y_Idx + N_box /2))
            CTT_Max = maxval(ctt(x_Idx_1:x_Idx_2,y_Idx_1:y_Idx_2))


            rain_column ( x_idx,y_idx) = min(dH_Max,(CTT_Max - CTT( x_idx,y_idx)) / Lapse_Rate + dH_0)


         end do x_loop
      end do y_loop

   end subroutine compute_rain_column

end module dncomp_precip_mod
