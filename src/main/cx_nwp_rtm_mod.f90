module cx_nwp_rtm_mod
  use RTM_COMMON_MOD, only: &
    NLevels_Rtm
    
  use univ_fp_comparison_mod, only: operator(.NEfp.)
  
  use CX_SCIENCE_TOOLS_MOD , only: &
       VAPOR  
    
    use NWP_COMMON_MOD, only: &
      NWP &
      , delta_t_inversion &
      , p_inversion_min 
   
   use CONSTANTS_MOD, only: &
      Real4 &
      , Int4 &
      , Int1 
      
    use PIXEL_COMMON_MOD, only: &
        NWP_PIX  
   
    use NUMERICAL_ROUTINES_MOD , only: &
       LOCATE
    implicit none
    
    private
    
    public:: convert_atmos_prof_nwp_rtm
    public:: map_nwp_rtm
    public:: create_temp_nwp_vectors
    public:: destroy_temp_nwp_vectors
    public :: Wvmr_Nwp
    
 integer, dimension(:), allocatable, save:: k_Nwp_Rtm
    real, dimension(:), allocatable, save:: x_Nwp_Rtm 
  integer, dimension(NLevels_Rtm) :: k_Rtm_Nwp
  real, dimension(NLevels_Rtm) :: x_Rtm_Nwp 
  
  
  real, dimension(:), allocatable,  save :: Wvmr_Nwp
  
contains

   !====================================================================
   ! subroutine Name: CREATE_TEMP_NWP_VECTORS
   !
   ! Function:
   !   create and initialize NWP vectors used for RTM calcs
   !
   !====================================================================
   subroutine CREATE_TEMP_NWP_VECTORS()

      allocate( Wvmr_Nwp(NWP%Nlevels), &
                k_Nwp_Rtm(NWP%Nlevels), &
                x_Nwp_Rtm(NWP%Nlevels))

      Wvmr_Nwp = 0.0
      k_Nwp_Rtm = 0.0
      x_Nwp_Rtm = 0.0

   end subroutine CREATE_TEMP_NWP_VECTORS
   
   
   !====================================================================
   ! subroutine Name: DESTROY_TEMP_NWP_VECTORS
   !
   ! Function:
   !   destroy and initialize NWP vectors used for RTM calcs
   !
   !====================================================================
   subroutine DESTROY_TEMP_NWP_VECTORS()

      deallocate(Wvmr_Nwp,k_Nwp_Rtm, x_Nwp_Rtm)

   end subroutine DESTROY_TEMP_NWP_VECTORS

   !====================================================================
   ! subroutine Name: MAP_NWP_RTM
   !
   ! Function:
   !   develops the mapping of the NWP and RTM profile Levels
   !
   !====================================================================
   subroutine MAP_NWP_RTM(Num_Levels_Nwp_Profile, &
                        Press_Profile_Nwp, &
                        Num_Levels_Rtm_Profile, &
                        Press_Profile_Rtm)

      integer, intent(in):: Num_Levels_Nwp_Profile
      integer, intent(in):: Num_Levels_Rtm_Profile
      real, intent(in), dimension(:)::  Press_Profile_Nwp
      real, intent(in), dimension(:)::  Press_Profile_Rtm

      integer:: k
      integer:: k_temp 
    
     
      !--- map NWP Levels to RTM Levels
      do k = 1, Num_Levels_Nwp_Profile
     
         !--- locate the nwp Level within the Rtm Levels 
         !--- P_Std_Nwp(k) should fall between Rtm Levels k_temp and k_temp +1
         call LOCATE(Press_Profile_Rtm, Num_Levels_Rtm_Profile, Press_Profile_Nwp(k), k_temp)
        
         k_Nwp_Rtm(k) = min(Num_Levels_Rtm_Profile-1,max(1,k_temp))

         !-- compute linear weighting factor
         x_Nwp_Rtm(k) = (Press_Profile_Nwp(k) - Press_Profile_Rtm(k_Nwp_Rtm(k))) / &
                      (Press_Profile_Rtm(k_Nwp_Rtm(k)+1) - Press_Profile_Rtm(k_Nwp_Rtm(k))) 
      end do

      !--- map RTM Levels to NWP Levels
      do k = 1, Num_Levels_Rtm_Profile
      
         !--- locate the Rtm Level within the Nwp Levels 
         !--- Press_Profile_Rtm(k) should fall between Nwp Levels k_temp and k_temp +1
         call LOCATE(Press_Profile_Nwp, Num_Levels_Nwp_Profile, Press_Profile_Rtm(k), k_temp)
         k_Rtm_Nwp(k) = min(Num_Levels_Nwp_Profile-1,max(1,k_temp))

         !-- compute linear weighting factor
         x_Rtm_Nwp(k) = (Press_Profile_Rtm(k) - Press_Profile_Nwp(k_Rtm_Nwp(k))) / &
                      (Press_Profile_Nwp(k_Rtm_Nwp(k)+1) - Press_Profile_Nwp(k_Rtm_Nwp(k)))
     
      end do 
      
   end subroutine MAP_NWP_RTM
   
   !====================================================================
   ! subroutine Name: CONVERT_ATMOS_PROF_NWP_RTM
   !
   ! Description:
   ! This routine interpolate the NWP profiles to profiles with the
   ! vertical spacing defined by the RTM model.  It operates on profiles
   ! stored in this module
   !
   ! INPUTS:
   !
   ! Highest_Level_Rtm_Nwp - highest Rtm Level that is below the highest nwp Level
   ! Lowest_Level_Rtm_Nwp - lowest Rtm Level that is above the lowest nwp Level
   ! Sfc_Level_Rtm - lowest Rtm Level above the surface
   ! P_near_Sfc_Nwp - lowest standard nwp Level above surface pressure
   !
   !====================================================================
   subroutine CONVERT_ATMOS_PROF_NWP_RTM(Num_Levels_NWP_Profile, &
                                       Surface_Level_Nwp, &
                                       Surface_Elevation_Nwp, &
                                       Air_Temperature_Nwp, &
                                       Surface_Rh_Nwp, &
                                       Surface_Pressure_Nwp, &
                                       Press_Profile_Nwp, &
                                       T_Profile_Nwp, &
                                       Z_Profile_Nwp, &
                                       Wvmr_Profile_Nwp, &
                                       Ozmr_Profile_Nwp, &
                                       Num_Levels_Rtm_Profile, &
                                       Press_Profile_Rtm, &
                                       T_Profile_Rtm, &
                                       Z_Profile_Rtm, &
                                       Wvmr_Profile_Rtm, &
                                       Ozmr_Profile_Rtm, &
                                       T_Std_Profile_Rtm, &
                                       Wvmr_Std_Profile_Rtm, &
                                       Ozmr_Std_Profile_Rtm)
                                        

      integer, intent(in):: Num_Levels_Nwp_Profile
      integer(kind=int1), intent(in):: Surface_Level_Nwp 
      real, intent(in):: Surface_Elevation_Nwp 
      real, intent(in):: Air_Temperature_Nwp 
      real, intent(in):: Surface_Rh_Nwp 
      real, intent(in):: Surface_Pressure_Nwp 
      real, intent(in), dimension(:):: Press_Profile_Nwp 
      real, intent(in), dimension(:):: T_Profile_Nwp 
      real, intent(in), dimension(:):: Z_Profile_Nwp 
      real, intent(in), dimension(:):: Wvmr_Profile_Nwp 
      real, intent(in), dimension(:):: Ozmr_Profile_Nwp 
      integer, intent(in):: Num_Levels_Rtm_Profile
      real, intent(in), dimension(:):: Press_Profile_Rtm
      real, intent(out), dimension(:):: T_Profile_Rtm
      real, intent(out), dimension(:):: Z_Profile_Rtm
      real, intent(out), dimension(:):: Wvmr_Profile_Rtm
      real, intent(out), dimension(:):: Ozmr_Profile_Rtm
      real, intent(in), dimension(:):: T_Std_Profile_Rtm
      real, intent(in), dimension(:):: Wvmr_Std_Profile_Rtm
      real, intent(in), dimension(:):: Ozmr_Std_Profile_Rtm

      integer:: k
      integer:: Lowest_Level_Rtm_Nwp
      integer:: Highest_Level_Rtm_Nwp
      integer:: Sfc_Level_Rtm
      real:: dT_dP_near_Sfc
      real:: dWvmr_dP_near_Sfc
      real:: dZ_dP_near_Sfc
      real:: P_near_Sfc_Nwp
      real:: Wvmr_Sfc
      real:: es
      real:: e
      real:: T_Offset

      !--- initialize indices
      Lowest_Level_Rtm_Nwp = Num_Levels_Rtm_Profile
      Highest_Level_Rtm_Nwp = 1
      Sfc_Level_Rtm = Num_Levels_Rtm_Profile
      P_Near_Sfc_Nwp = Press_Profile_Nwp(Surface_Level_Nwp)
   
      !--- make Wvmr at sfc
      es = VAPOR(Air_Temperature_Nwp)
      e = Surface_Rh_Nwp * es / 100.0
      Wvmr_Sfc = 1000.0*0.622 * (e / (Surface_Pressure_Nwp - e))  !(g/kg)

      !--- determine some critical Levels in the Rtm profile
      do k = 1, Num_Levels_Rtm_Profile
         if (Press_Profile_Rtm(k) > Press_Profile_Nwp(1)) then
            Highest_Level_Rtm_Nwp = k
            exit
         endif
      enddo

      do k = Num_Levels_Rtm_Profile,1,-1
         if (Press_Profile_Rtm(k) < Surface_Pressure_Nwp) then
            Sfc_Level_Rtm = k
            exit
         endif
      enddo

      do k = Num_Levels_Rtm_Profile,1,-1
         if (Press_Profile_Rtm(k) < P_Near_Sfc_Nwp) then
            Lowest_Level_Rtm_Nwp = k
            exit
         endif
      enddo

      !--- compute T and Wvmr lapse rate near the surface
      dT_dP_near_Sfc = 0.0
      dZ_dP_near_Sfc = 0.0
      dWvmr_dP_near_Sfc = 0.0

      if (Surface_Pressure_Nwp .NEfp. Press_Profile_Nwp(Num_Levels_Nwp_Profile)) then
         dT_dP_near_Sfc =  &
             (T_Profile_Nwp(Surface_Level_Nwp) - Air_Temperature_Nwp)/ &
             (Press_Profile_Nwp(Surface_Level_Nwp) - Surface_Pressure_Nwp)
         dWvmr_dP_near_Sfc =  &
             (Wvmr_Profile_Nwp(Surface_Level_Nwp) - Wvmr_Sfc)/ &
             (Press_Profile_Nwp(Surface_Level_Nwp) - Surface_Pressure_Nwp)
      else
         dT_dP_near_Sfc =  &
             (T_Profile_Nwp(Surface_Level_Nwp-1) - Air_Temperature_Nwp)/ &
             (Press_Profile_Nwp(Surface_Level_Nwp-1) - Surface_Pressure_Nwp)
         dWvmr_dP_near_Sfc =  &
             (Wvmr_Profile_Nwp(Surface_Level_Nwp-1) - Wvmr_Sfc)/ &
             (Press_Profile_Nwp(Surface_Level_Nwp-1) - Surface_Pressure_Nwp)
      end if
   
      if (Press_Profile_Nwp(Surface_Level_Nwp-1) .NEfp. Press_Profile_Nwp(Surface_Level_Nwp)) then
         dZ_dP_near_Sfc =  &
             (Z_Profile_Nwp(Surface_Level_Nwp-1) - Z_Profile_Nwp(Surface_Level_Nwp))/ &
             (Press_Profile_Nwp(Surface_Level_Nwp-1) - Press_Profile_Nwp(Surface_Level_Nwp))
       
      end if

      ! --- compute temperature offset between standard and NWP profiles at top
      !--- this will be added to the standard profiles
      T_Offset = T_Profile_Nwp(1) - T_Std_Profile_Rtm(Highest_Level_Rtm_Nwp)

      !--- for Rtm Levels above the highest nwp Level, use Rtm standard
      !values
      do k = 1,Highest_Level_Rtm_Nwp-1
          Z_Profile_Rtm(k) = Z_Profile_Nwp(1)
          T_Profile_Rtm(k) = T_Std_Profile_Rtm(k) + T_Offset
          Wvmr_Profile_Rtm(k) = Wvmr_Std_Profile_Rtm(k)
          Ozmr_Profile_Rtm(k) = Ozmr_Std_Profile_Rtm(k)
      end do

      !--- Rtm Levels within standard nwp Levels above the surface
      do k = Highest_Level_Rtm_Nwp, Lowest_Level_Rtm_Nwp
          T_Profile_Rtm(k) = T_Profile_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (T_Profile_Nwp(k_Rtm_Nwp(k)+1) - T_Profile_Nwp(k_Rtm_Nwp(k)))
          Z_Profile_Rtm(k) = Z_Profile_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (Z_Profile_Nwp(k_Rtm_Nwp(k)+1) - Z_Profile_Nwp(k_Rtm_Nwp(k)))
          Wvmr_Profile_Rtm(k) = Wvmr_Profile_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (Wvmr_Profile_Nwp(k_Rtm_Nwp(k)+1) - Wvmr_Profile_Nwp(k_Rtm_Nwp(k)))

          Ozmr_Profile_Rtm(k) = Ozmr_Profile_Nwp(k_Rtm_Nwp(k)) + x_Rtm_Nwp(k) *  &
                     (Ozmr_Profile_Nwp(k_Rtm_Nwp(k)+1) - Ozmr_Profile_Nwp(k_Rtm_Nwp(k)))
      end do

      !--- Rtm Levels that are below the lowest nwp Level but above the surface
      do k = Lowest_Level_Rtm_Nwp+1, Sfc_Level_Rtm
           T_Profile_Rtm(k) = Air_Temperature_Nwp + dT_dP_near_Sfc * &
                      (Press_Profile_Rtm(k) - Surface_Pressure_Nwp)
           Wvmr_Profile_Rtm(k) = Wvmr_Sfc + dWvmr_dP_near_Sfc * &
                      (Press_Profile_Rtm(k) - Surface_Pressure_Nwp)
           Z_Profile_Rtm(k) = Surface_Elevation_Nwp + dZ_dP_near_Sfc * &
                      (Press_Profile_Rtm(k) - Surface_Pressure_Nwp)
           Ozmr_Profile_Rtm(k) = Ozmr_Profile_Nwp(Num_Levels_Nwp_Profile)
      end do

      !--- Rtm Levels below the surface
      do k = Sfc_Level_Rtm +1, Num_Levels_Rtm_Profile
         T_Profile_Rtm(k) = Air_Temperature_Nwp
         Z_Profile_Rtm(k) = Surface_Elevation_Nwp + &
                            dZ_dP_near_Sfc * (Press_Profile_Rtm(k) - Surface_Pressure_Nwp)
         Wvmr_Profile_Rtm(k) = Wvmr_Sfc
         Ozmr_Profile_Rtm(k) = Ozmr_Std_Profile_Rtm(k)
      end do

      !--- if using NCEP reanalysis which has no ozone profile, use default
      if (NWP_PIX%Nwp_Opt == 2) then
         Ozmr_Profile_Rtm = Ozmr_Std_Profile_Rtm
      end if

   end subroutine CONVERT_ATMOS_PROF_NWP_RTM

end module cx_nwp_rtm_mod
