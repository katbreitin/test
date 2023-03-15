!  $Header: https://svn.ssec.wisc.edu/repos/cloud_team_clavrx/trunk/main_src/dnb_retrievals_mod.f90 3082 2018-12-17 17:53:19Z mfoster $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: dnb_retrievals_mod.f90 (src)
!       dnb_retrievals_mod (program)
!
! PURPOSE: this program computes lunar eflectance
!
! DESCRIPTION:
!
! AUTHORS:
!  Steven Miller , CIRA
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
!  Denis Botambekov, CIMSS, denis.botambekov@ssec.wisc.edu
!  William Straka, CIMSS, wstraka@ssec.wisc.edu
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
!HISTORY:
!   02/20/2015: updated coeffcients sent from Steve on 26 Jan 2015 (AW)
!
!--------------------------------------------------------------------------------------
module dnb_retrievals_mod
   use FILE_UTILS, only: file_test
   use univ_kind_defs_mod, only: f4, i2_B, in, f8, i4_B

   use univ_fp_comparison_mod, only: operator(.EQfp.), operator(.GEfp.)

   implicit none
   private
   public :: compute_lunar_reflectance
   public :: lunar_reflectance_nasa_data_adjustment
   private:: index_in_vector

   integer, parameter, private :: INT2 = selected_int_kind(3)
   integer, parameter, private :: INT4 = selected_int_kind(8)
   integer, parameter, private :: REAL4 = selected_real_kind(6,37)

   ! - paramaters
   ! - physical
   real(f8) , parameter, private :: MEAN_EARTH_SUN_DIST = 149598022.6071
   real(f8) , parameter, private :: MEAN_EARTH_MOON_DIST = 384400.0
   real , parameter , private :: EARTH_RADIUS_KM = 6378.140
   real , parameter , private:: SRF_INTEG = 0.32560294 ! integral of the DNB sensor response function (micron)
   real , parameter , private:: ASTRO_DARK_THRESH = 109.0 ! Sun 19 degrees or more below horizon
   real , parameter , private :: MIN_LUNAR_IRRAD_DNB = 1.0e-5 ! W/m^2 = 1.0e-09 W/cm^2, threshold for doing calcs
   !- other params
   real , parameter, private :: PI = ACOS(-1.)
   real , parameter, private ::  DTOR = PI / 180.
   integer, parameter, private :: FLOAT_SIZE  = 4
   integer, parameter, private :: DOUBLE_SIZE = 8

contains

   subroutine compute_lunar_reflectance ( &
               rad_chdnb_input &
             & , solzen &
             & , lunzen &
             & , start_year &
             & , month &
             & , day_of_month &
             & , start_time &
             & , lunar_phase_angle_topo &
             & , ancil_data_dir &
             & , ref_chdnb_lunar)


      implicit none

      ! - input
      real , intent(in) :: rad_chdnb_input(:,:)
      real , intent(in) :: solzen(:,:)
      real , intent(in) :: lunzen(:,:)
      integer(i4_B), intent(in) :: start_time
      integer(i4_B), intent(in) :: start_year
      integer(i4_B), intent(in) :: month
      integer(i4_B), intent(in) :: day_of_month
      real(f8) , intent(in) :: lunar_phase_angle_topo
      character ( len = * ), intent(in) :: ancil_data_dir

      ! - output
      !---akh-- real, intent ( out ) ,allocatable :: ref_chdnb_lunar(:,:)
      real, intent ( out ) :: ref_chdnb_lunar(:,:)

      !********
      ! NEW 3/26/2014 based on Gauss curve-fits to Obs/Mod ratio data between -120(wax) and 120(wane) degrees
      ! Note: "lpds" means "lunar phase, degrees, signed"
      real(f8)::  lpds

      ! USE 4/10/2015 more precision

      ! WAX VIS08

      real(f8):: waxp1=5.2228373e-12
      real(f8):: waxp2=2.2410515e-9
      real(f8):: waxp3=3.7994231e-7
      real(f8):: waxp4=3.1637454e-5
      real(f8):: waxp5=1.3079265e-3
      real(f8):: waxp6=2.3328238e-2
      real(f8):: waxp7=1.1448359


      ! WANE VIS08

      real(f8):: wanp1=6.3455594e-12
      real(f8):: wanp2=-2.6095939e-9
      real(f8):: wanp3=4.2557303e-7
      real(f8):: wanp4=-3.4087847e-5
      real(f8):: wanp5=1.3562948e-3
      real(f8):: wanp6= -2.5037150e-2
      real(f8):: wanp7=1.1450824


      real, dimension(:,:) , allocatable :: rad_chdnb

      real(f8) :: yyyymmddhh
      real(f8) :: lunar_phase_angle_geo
      real(f8) :: cos_phase_angle
      real(f8) :: curr_earth_sun_dist
      real(f8) :: curr_earth_moon_dist
      real(f8) :: curr_mean_irrad
      real(f8) :: cos_weighted_irrad


      integer( kind = 4 ) , parameter :: num_irrad_tabvals = 181
      integer( kind = 4 ) , parameter :: num_dist_tabvals = 184080
      real, allocatable :: lunar_irrad_lut (:,:)
      real(f8), allocatable :: dist_phase_lut(:,:)
      real(f8), allocatable :: phase_array(:)


      real :: hour_fraction
      real(f8) :: phase_fraction
      real :: minute
      real :: hour
      real(f8) :: lunar_irrad_dnb
      integer :: i
      integer :: j
      integer :: dtg_index
      integer ::  irrad_index
      double precision :: denorm1
      double precision :: denorm2
      double precision :: denorm3
      double precision :: denorm_factor

      logical :: dnb_verbose = .false.
      !logical :: dnb_verbose = .true.
      character(len=1020) :: lunar_irrad_file
      character(len=1020) :: distance_table_file
      integer :: num_pix , num_elem

      logical :: is_waning
      real(f8) :: phase_albedo_correction_factor = 1.



      ! --- executable --------------------------------------

      lunar_irrad_file     = trim(ancil_data_dir)//'static/dnb_ancils/lunar_irrad_Mean_DNB.bin'
      if ( .not. file_test ( trim(lunar_irrad_file) ) ) then
         print* , 'lunar irradiance file missing ', lunar_irrad_file
         return
      end if

      distance_table_file  = trim(ancil_data_dir)//'static/dnb_ancils/DIST_2010-2030_double.bin'
      if ( .not. file_test ( trim(distance_table_file) ) ) then
         print* , 'lunar irradiance file missing ', distance_table_file
         return
      end if


      num_pix = ubound(rad_chdnb_input,1)
      num_elem = ubound(rad_chdnb_input,2)

      allocate ( lunar_irrad_LUT(2,num_irrad_tabvals) )
      allocate ( dist_phase_LUT(4,num_dist_tabvals))
      allocate ( phase_array(num_irrad_tabvals) )
!---  allocate ( ref_chdnb_lunar(num_pix,num_elem) )
      allocate ( rad_chdnb(num_pix, num_elem))

      ! - read LUT values   for Sun Earth Moon geometry
      open (unit=1,file=trim(lunar_irrad_file),status="old",action="read",&
            access="direct",form="unformatted",recl=float_size*2*num_irrad_tabvals)
      read (unit=1,rec=1) lunar_irrad_lut
      close (1)

      open (unit=1,file=trim(distance_table_file),status="old",action="read",&
             access="direct",form="unformatted",recl=double_size*4*num_dist_tabvals)
      read (unit=1,rec=1) dist_phase_lut
      close (1)


      ! -  3. compute toa downwelling lunar irradiance (lunar_irrad_dnb) for current date/time
      minute = mod(start_time/1000./60., 60.)
      hour_fraction = minute / 60.0
      hour = floor(start_time/1000./60./60.)
      yyyymmddhh = start_year * 1000000 + month*10000+ day_of_month *100 + int(hour)

      dtg_index = index_in_vector(dist_phase_lut(1,:),num_dist_tabvals,yyyymmddhh)



      dtg_index = min(num_dist_tabvals,dtg_index)
      dtg_index = max(1,dtg_index)

      lunar_phase_angle_geo    = dist_phase_lut(2,dtg_index)  &
                        & +  hour_fraction &
                        & * (dist_phase_lut(2,dtg_index + 1) &
                        & - dist_phase_lut(2,dtg_index))

      curr_earth_sun_dist  = dist_phase_lut(3,dtg_index)  &
                        & + hour_fraction &
                        & * (dist_phase_lut(3,dtg_index + 1) &
                        & - dist_phase_lut(3,dtg_index))

      curr_earth_moon_dist = dist_phase_lut(4,dtg_index)  &
                        & + hour_fraction &
                        & * (dist_phase_lut(4,dtg_index + 1) &
                        & - dist_phase_lut(4,dtg_index))

      if (dnb_verbose) then

         print *,''
         print *,'compare (these two values should be about the same):'
         print *, 'lunar_phase (from viirs granule) = ',lunar_phase_angle_topo
         print *, 'lunar_phase_angle_geo (from lut) = ',lunar_phase_angle_geo
         print *,''
         print *,''
         print *,'dist_phase_lut(:,dtg_index) = ',dist_phase_lut(:,dtg_index)
         print *,'dist_phase_lut(:,dtg_index+1) = ',dist_phase_lut(:,dtg_index+1)
         print *,'hour_fraction = ',hour_fraction
         print *,'lunar_phase_angle_geo = ',lunar_phase_angle_geo
         print *,'curr_earthsun_dist = ',curr_earth_sun_dist
         print *,'curr_earthmoon_dist = ',curr_earth_moon_dist
         print *,''
         print *,'solar angle min max ', minval(solzen),maxval(solzen)
      end if

      ! b) interpolate lunar_irrad_lut() to get current mean-geometry lunar irradiance pre-convolved to dnb srf
      !   use topo phase angle
      phase_fraction  = lunar_phase_angle_topo - int(lunar_phase_angle_topo)
      phase_array     = lunar_irrad_lut(1,:)
      if (lunar_phase_angle_topo < phase_array(1) .or.  &
           lunar_phase_angle_topo > phase_array(num_irrad_tabvals)) then
         ! lunar_phase_angle_topo is out of bounds; set canary value and return
         ref_chdnb_lunar = -999.0
         return
      end if
      irrad_index     = index_in_vector(phase_array,num_irrad_tabvals,lunar_phase_angle_topo)
      curr_mean_irrad = lunar_irrad_lut(2,irrad_index) + &
                     &   phase_fraction*(lunar_irrad_lut(2,irrad_index+1)-lunar_irrad_lut(2,irrad_index))

      deallocate ( phase_array)
      if (dnb_verbose) then
         print *,'phase_fraction = ',phase_fraction
         print *,'irrad_index = ',irrad_index
         print *,'lunar_irrad_lut(:,irrad_index) = ',lunar_irrad_lut(:,irrad_index)
         print *,'lunar_irrad_lut(:,irrad_index+1) = ',lunar_irrad_lut(:,irrad_index+1)
         print *,'curr_mean_irrad = ',curr_mean_irrad
      end if

      deallocate ( lunar_irrad_lut )

      !  c) define denormalization parameters to scale irradiance to current sun/earth/moon geometry
      !   use geo phase angle
      cos_phase_angle = cos(lunar_phase_angle_geo * DTOR)
      denorm1 = mean_earth_sun_dist ** 2 + mean_earth_moon_dist ** 2 + &
           2.0 * mean_earth_moon_dist * mean_earth_sun_dist * cos_phase_angle
      denorm2 = curr_earth_sun_dist ** 2 + curr_earth_moon_dist ** 2 + &
           2.0 * curr_earth_moon_dist * curr_earth_sun_dist * cos_phase_angle
      denorm3 = ((mean_earth_moon_dist - earth_radius_km ) / (curr_earth_moon_dist - earth_radius_km)) ** 2.0
      denorm_factor = (denorm1 / denorm2 ) * denorm3


      !  d) denormalize mean-geometry irradiance to current geometry
      !     also, convert from mw/m^2-um to w/m^2 (divide by 1000 and multiply by srf_integ)
      lunar_irrad_dnb = curr_mean_irrad * denorm_factor * (SRF_INTEG * 1.0e-03)

      if (dnb_verbose) then
         print *,'curr_mean_irrad (mw/m^2-micron)= ',curr_mean_irrad
         print *,'denorm_factor = ',denorm_factor
         print *,'--> dnb band-integrated lunar irradiance (w/m^2)= ',lunar_irrad_dnb
      end if


      !********
      ! PHASE ANGLE BIAS CORRECTION
      !********
      !  e) Compute phase-angle-dependent correction term to account for albedo variation
      !     Based on curve fits to the Obs/Modeled ratios
      !     Make sure that we are using the appropriate set of coefficients
      !     Expansion works on the signed lunar phase angle values (degrees)


      ! 1 determine waxing or waning and adjust lunar_phase

      lpds = lunar_phase_angle_topo
      is_waning = dist_phase_LUT(2,dtg_index ) > dist_phase_LUT(2,dtg_index + 1)
      if ( is_waning ) lpds = lpds * ( -1 )

      deallocate ( dist_phase_lut )

      if ( abs ( lpds ) > 120.0 ) then
         phase_albedo_correction_factor = 1.
      else if ( lpds < 0 ) then
         phase_albedo_correction_factor = waxp1*lpds**6 + waxp2*lpds**5 + waxp3*lpds**4 + &
                                          waxp4*lpds**3  +waxp5*lpds**2 + waxp6*lpds + waxp7
      else
         phase_albedo_correction_factor = wanp1*lpds**6 + wanp2*lpds**5 + wanp3*lpds**4 + &
                                          wanp4*lpds**3  +wanp5*lpds**2 + wanp6*lpds + wanp7
      end if

      lunar_irrad_dnb = lunar_irrad_dnb * phase_albedo_correction_factor

      rad_chdnb= rad_chdnb_input  * 1.0e+04
      ref_chdnb_lunar = -999.0
      do i=1,num_pix
         do j=1,num_elem

            if (rad_chdnb( i , j ) < 0 .or. &
                  & solzen( i , j ) < astro_dark_thresh .or. &
                  & lunzen( i , j ) > 90.0) cycle

            if (rad_chdnb ( i , j ) > -1.0) ref_chdnb_lunar( i , j ) = 0.0

            cos_weighted_irrad = cos(lunzen( i , j ) * DTOR) * lunar_irrad_dnb

            if (cos_weighted_irrad > MIN_LUNAR_IRRAD_DNB) then
               ref_chdnb_lunar( i , j ) = 100.0 * ( PI * rad_chdnb( i , j ) ) &
                 / ( real(cos_weighted_irrad,f4) )

               if (ref_chdnb_lunar( i , j ) > 150.0) then
                  ref_chdnb_lunar( i , j ) = 150.0
               end if

            else
               ref_chdnb_lunar( i , j ) = 0.0
            end if

         end do !j
      end do !i

      deallocate ( rad_chdnb)

   end subroutine compute_lunar_reflectance

   !
   !
   !
    subroutine lunar_reflectance_nasa_data_adjustment (rfl_lunar)
      real , intent(inout) :: rfl_lunar(:,:)
      real, parameter :: Dnb_Coef(3) = [-0.118767,0.962452,-0.000144502]

      where( rfl_lunar .GEfp. 0.00)
         rfl_lunar = Dnb_Coef(1) + &
                            Dnb_Coef(2) * rfl_lunar + &
                            Dnb_Coef(3) * rfl_lunar **2
      endwhere


    end subroutine lunar_reflectance_nasa_data_adjustment

   !----------------------------------------------------------------------
   !
   !
   !----------------------------------------------------------------------
   integer function index_in_vector(xx,n,x)
         integer, intent(in)::n
         real(f8), intent(in)::x
         real(f8),dimension(:),intent(in)::xx
         integer:: i,jl,jm,ju
         jl=0
         ju=n+1

         do i=1,2*n
            if (ju-jl <= 1) then
               exit
            end if
            jm=(ju+jl)/2
            if ( (xx(n) >= xx(1)) .eqv. (x >= xx(jm)) ) then
               jl=jm
            else
               ju=jm
            end if
         end do

         if (x .EQfp. xx(1)) then
            index_in_vector =1
         else if (x .EQfp. xx(n)) then
            index_in_vector = n-1
         else
            index_in_vector = jl
         end if

   end function index_in_vector




end module dnb_retrievals_mod
