module cx_atmospheric_correction_vis_mod

 use PIXEL_COMMON_MOD
 use univ_fp_comparison_mod, only: operator(.EQfp.), operator(.NEfp.)
  use CLAVRX_MESSAGE_MOD, only: MESG, verb_lev
  use SURFACE_PROPERTIES_MOD, only: &
        ch1_sfc_alb_umd &
      , ch2_sfc_alb_umd &
      , ch6_sfc_alb_umd &
      , ch1_snow_sfc_alb_umd &
      , ch2_snow_sfc_alb_umd &
      , ch5_snow_sfc_alb_umd &
      , ch6_snow_sfc_alb_umd &
      , ch7_snow_sfc_alb_umd 
 
 
 use CONSTANTS_MOD
 
 public::  &
            SETUP_SOLAR_RTM, &
            COMPUTE_CLEAR_SKY_SCATTER, &
            ATMOS_CORR
       
       
       

contains


   !==============================================================================
   !
   !==============================================================================
   subroutine SETUP_SOLAR_RTM ( WMO_Id )

      integer, intent(in):: WMO_Id
      real, parameter:: alpha = 0.5
      real, parameter:: Tau_Aer_1um = 0.0 !0.10

      !------------------------------------------------------------------------
      ! initialize gas coefficients with default values
      !------------------------------------------------------------------------
      include 'default_solar_terms.inc'

      !---- define channel aerosol properties
      Solar_Rtm%Tau_Aer = 0.0
      Solar_Rtm%Wo_Aer = 0.8
      Solar_Rtm%G_Aer = 0.6

      ! tau(lambda) = tau(lambda=1) * lambda**(-alpha)
      Solar_Rtm%Tau_Aer(1) = Tau_Aer_1um * (0.65**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(2) = Tau_Aer_1um * (0.86**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(3) = Tau_Aer_1um * (0.47**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(5) = Tau_Aer_1um * (0.55**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(6) = Tau_Aer_1um * (1.60**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(7) = Tau_Aer_1um * (2.15**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(8) = Tau_Aer_1um * (0.41**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(9) = Tau_Aer_1um * (0.44**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(10) = Tau_Aer_1um * (0.49**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(11) = Tau_Aer_1um * (0.53**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(12) = Tau_Aer_1um * (0.55**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(13) = Tau_Aer_1um * (0.67**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(14) = Tau_Aer_1um * (0.68**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(15) = Tau_Aer_1um * (0.75**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(16) = Tau_Aer_1um * (0.87**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(17) = Tau_Aer_1um * (0.90**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(18) = Tau_Aer_1um * (0.93**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(19) = Tau_Aer_1um * (0.94**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(20) = Tau_Aer_1um * (3.75**(-1.0*alpha))
      Solar_Rtm%Tau_Aer(26) = Tau_Aer_1um * (1.38**(-1.0*alpha))


      select case(WMO_Id)

      include 'avhrr_solar_terms.inc'
      include 'goes_ip_solar_terms.inc'
      include 'viirs_solar_terms.inc'
      include 'seviri_solar_terms.inc'
      include 'modis_solar_terms.inc'
      include 'abi_solar_terms.inc'
      include 'ahi_solar_terms.inc'
      include 'coms_solar_terms.inc'
      include 'mtsat_solar_terms.inc'
      include 'metimage_solar_terms.inc'

!     case(174)     ! AHI-9 VALUES NEED TO BE UPDATED

      case default
         call MESG("WARNING: Solar transmission terms are not available, using default", level = verb_lev % DEFAULT)
      end select

   end subroutine SETUP_SOLAR_RTM
!------------------------------------------------------------
! atmospheric correction for reflectance channels
!
! note, original optical depths taken from 
!       Ignatov and Stowe, 2002, vol 59, JAS, pg 313-334
!
! AKH nows runs modtran to generate these numbers
!-----------------------------------------------------------
subroutine ATMOS_CORR(Line_Idx_Min,Num_Lines)

   integer, intent(in):: Line_Idx_Min
   integer, intent(in):: Num_Lines
   integer:: Line_Idx
   integer:: Elem_Idx
   integer:: Line_Idx_Max
   integer:: Elem_Idx_Max
   integer:: Elem_Idx_Min
   integer:: Num_Elements
   integer:: Chan_Idx

   real:: Ref_ss
   real:: Albedo_View
   real:: Albedo_Sun
   real:: Trans_Total
   real:: Tau_Total
   real:: Tau_Ray
   real:: Tau_Gas
   real:: Tau_H2O
   real:: Tau_Aer
   real:: Wo_Aer
   real:: G_Aer
   real:: Tpw_Ac
   real:: Zc, Pc

   real:: Source_Zen
   real:: Cos_Source_Zen
   real:: Airmass_Factor
   real:: Scattering_Angle

   !--- compute gas terms
   real, parameter :: h_h2o = 2000.0 !m
   real, parameter :: P_Sfc = 1013.25

   Elem_Idx_Min = 1
   Num_Elements = Image%Number_Of_Elements
   Elem_Idx_Max = Elem_Idx_Min + Num_Elements - 1
   Line_Idx_Max = Line_Idx_Min + Num_Lines - 1

  !---------- loop over scan lines
  line_loop: do Line_Idx = Line_Idx_Min, Line_Idx_Max

  element_loop: do Elem_Idx = Elem_Idx_Min, Elem_Idx_Max

     !--- check for bad individual pixels
     if (Bad_Pixel_Mask(Elem_Idx,Line_Idx) == sym%YES) cycle

     channel_loop: do Chan_Idx = 1,Nchan_Clavrx

       if (Sensor%Chan_On_Flag_Default(Chan_Idx) == sym%NO) cycle

       if (Chan_Idx /= 1 .and. Chan_Idx /= 2 .and. Chan_Idx /= 3 .and. &
           Chan_Idx /= 4 .and. Chan_Idx /= 5 .and. Chan_Idx /= 6 .and. &
           Chan_Idx /= 7 .and. Chan_Idx /= 26 .and. Chan_Idx /= 44) cycle

       !--- check for valid data
       if (Chan_Idx /= 44 .and. (ch(Chan_Idx)%Ref_Toa(Elem_Idx,Line_Idx) .EQfp. Missing_Value_Real4)) cycle
       if (Chan_Idx == 44) then
         if (ch(Chan_Idx)%Ref_Lunar_Toa(Elem_Idx,Line_Idx) .EQfp. Missing_Value_Real4) cycle
       endif

       !--- set source angle
       Source_Zen = Geo%SolZen(Elem_Idx,Line_Idx)
       Scattering_Angle = Geo%Scatangle(Elem_Idx,Line_Idx)
       if (Chan_Idx == 44) then
            Source_Zen = Geo%LunZen(Elem_Idx,Line_Idx)
            Scattering_Angle = Geo%Scatangle_Lunar(Elem_Idx,Line_Idx)
       endif

       !--- check for appropriate illumination
       if (Source_Zen >= 90.0) cycle

       !--- compute gas terms
       Tpw_Ac = NWP_PIX%Tpw(Elem_Idx,Line_Idx)
       if (Zc_Opaque_Cloud(Elem_Idx,Line_Idx) .NEfp. MISSING_VALUE_REAL4) then
          Zc = max(0.0,Zc_Opaque_Cloud(Elem_Idx,Line_Idx))
          Tpw_Ac = Tpw_Ac * exp(-Zc/h_h2o)
       endif

       Tau_H2O = Solar_Rtm%Tau_H2O_Coef(Chan_Idx,1) + Solar_Rtm%Tau_H2O_Coef(Chan_Idx,2)*Tpw_Ac +  &
                 Solar_Rtm%Tau_H2O_Coef(Chan_Idx,3)*(Tpw_Ac**2)

       Tau_Gas = max(0.0,Tau_H2O) + Solar_Rtm%Tau_O3(Chan_Idx) + Solar_Rtm%Tau_O2(Chan_Idx) &
               + Solar_Rtm%Tau_CO2(Chan_Idx) + Solar_Rtm%Tau_CH4(Chan_Idx)

       Tau_Ray = Solar_Rtm%Tau_Ray(Chan_Idx)
       Tau_Aer = Solar_Rtm%Tau_Aer(Chan_Idx)
       if (Pc_Opaque_Cloud(Elem_Idx,Line_Idx) .NEfp. MISSING_VALUE_REAL4) then
          Pc = min(1000.0,max(0.0,Pc_Opaque_Cloud(Elem_Idx,Line_Idx)))
          Tau_Ray = Tau_Ray * Pc_Opaque_Cloud(Elem_Idx,Line_Idx) / P_Sfc
          Tau_Aer = Tau_Aer * (Pc_Opaque_Cloud(Elem_Idx,Line_Idx) / P_Sfc)**2
       endif

       Wo_Aer = Solar_Rtm%Wo_Aer(Chan_Idx)
       G_Aer = Solar_Rtm%G_Aer(Chan_Idx)

       !------------------------------------------------------------------------
       ! select gas and surface reflectance parameters
       !------------------------------------------------------------------------
       select case (Chan_Idx)

       case(1)
         if (ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) .NEfp. Missing_Value_Real4) then
              Albedo_View = ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else
              Albedo_View = Ch1_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
              Albedo_View = Ch1_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(2)
         if (ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) .NEfp. Missing_Value_Real4) then
              Albedo_View = ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else
              Albedo_View = Ch2_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
              Albedo_View = Ch2_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(3)
         if (ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) .NEfp. Missing_Value_Real4) then
              Albedo_View = ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else
              Albedo_View = Ch1_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
              Albedo_View = Ch1_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(4)
         if (ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) .NEfp. Missing_Value_Real4) then
              Albedo_View = ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else
              Albedo_View = Ch1_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
              Albedo_View = Ch1_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(5)
         if (ch(5)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) .NEfp. Missing_Value_Real4) then
              Albedo_View = ch(5)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else
              Albedo_View = Ch6_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0 !Note there is no Ch5_Sfc_Alb_Umd
         endif
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
              Albedo_View = Ch5_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(6)
         if (ch(6)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) .NEfp. Missing_Value_Real4) then
              Albedo_View = ch(6)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else
              Albedo_View = Ch6_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
              Albedo_View = Ch6_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(7)
         if (ch(7)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) .NEfp. Missing_Value_Real4) then
              Albedo_View = ch(7)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else
              Albedo_View = Ch6_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0   !Note there is no Ch7_Sfc_Alb_Umd
         endif
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
              Albedo_View = Ch7_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(26)
         if (ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) .NEfp. Missing_Value_Real4) then
              Albedo_View = ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) / 100.0
         else
              Albedo_View = Ch1_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0   !Note there is no Ch26_Sfc_Alb_Umd
         endif
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
              Albedo_View = Ch1_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) / 100.0
         endif

       case(44)  !DNB - use mean of ch1 and ch2 for sfc reflectance
         if ((ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) .NEfp. Missing_Value_Real4) &
            .and. (ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx) .NEfp. Missing_Value_Real4)) then
              Albedo_View = 0.5*(ch(1)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx)+ch(2)%Sfc_Ref_White_Sky(Elem_Idx,Line_Idx)) / 100.0
         else
              Albedo_View = 0.5*(Ch1_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) + Ch2_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx))) / 100.0
         endif
         if (Sfc%Snow(Elem_Idx,Line_Idx) /= sym%NO_SNOW) then
              Albedo_View = 0.5*(Ch1_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx)) &
               + Ch2_Snow_Sfc_Alb_Umd(Sfc%Sfc_Type(Elem_Idx,Line_Idx))) / 100.0
         endif
       end select

       !-- assume lambertian
       Albedo_Sun = Albedo_View

       Cos_Source_Zen = Missing_Value_Real4
       if (Source_Zen >= 0.0 .and. Source_Zen < 90.0) Cos_Source_Zen = cos(Source_Zen*DTOR)

       Airmass_Factor = 1.0 / Geo%CosZen(Elem_Idx,Line_Idx) + &
                        1.0 / Cos_Source_Zen

       !--- compute atmospheric scattering
       call COMPUTE_CLEAR_SKY_SCATTER(Tau_Aer, &
                                     Wo_Aer, &
                                     G_Aer, &
                                     Tau_Ray, &
                                     Tau_Gas, &
                                     Scattering_Angle, &
                                     Geo%Coszen(Elem_Idx,Line_Idx), &
                                     Cos_Source_Zen, &
                                     Albedo_View, &
                                     Albedo_Sun, &
                                     Ref_ss)

       !--- compute total transmission for combining terms
       Tau_Total = Tau_Aer + Tau_Ray + Tau_Gas
       Trans_Total = exp(-Tau_Total*Airmass_Factor)

       if (Chan_Idx /= 44) then

         !--- compute atmospherically corrected reflectance (at sfc level)s
         ch(Chan_Idx)%Ref_Sfc(Elem_Idx,Line_Idx) = (ch(Chan_Idx)%Ref_Toa(Elem_Idx,Line_Idx) - Ref_ss) / Trans_Total

         !--- compute top of clear-sky atmosphere reflectance
         ch(Chan_Idx)%Ref_Toa_Clear(Elem_Idx,Line_Idx) = Ref_ss + Trans_Total*100.0*Albedo_View

       else

         !--- compute atmospherically corrected reflectance (at sfc level)s
         ch(Chan_Idx)%Ref_Lunar_Sfc(Elem_Idx,Line_Idx) = (ch(Chan_Idx)%Ref_Lunar_Toa(Elem_Idx,Line_Idx) - Ref_ss) / Trans_Total

         !--- compute top of clear-sky atmosphere reflectance
         ch(Chan_Idx)%Ref_Lunar_Toa_Clear(Elem_Idx,Line_Idx) = Ref_ss + Trans_Total*100.0*Albedo_View

       endif

    end do  channel_loop

   end do element_loop

 end do line_loop

end subroutine ATMOS_CORR


!------------------------------------------------------------
! compute the single scater and aerosol reflectance
! this assumes that the gas is mixed in with scattering
!
! Input: 
!  scatangle = scattering angle
!  coszen = cosine of viewing zenith angle
!  cossolzen = cosine of solar zenith angle
!  tau_ray = rayleigh optical depth
!  tau_gas = gaseous optical depth
!  tau_aer = aerosol optical depth
!  wo_aer = aerosol single scatter albedo
!  g_aer = aerosol asymmetry parameter
!  cloud_albedo_zen = cloud or surface albedo from illumination along coszen
!  cloud_albedo_solzen = cloud or surface albedo from illumination along
!  cosolzen
!  note- albedoes are (0-1)
!
! Output:
!   Ref_SS = single scatter reflectance from the layer above the 
!            surface or cloud at the top of atmosphere (0-110%)
!
! Reference:
! Correction of Rayleigh scattering effects in cloud optical thickness retrievals 
! Menghua Wang and Michael D. King 
! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 102, NO. D22, PAGES 25,915-25,926,
! NOVEMBER 27, 1997 
!-----------------------------------------------------------
subroutine COMPUTE_CLEAR_SKY_SCATTER(Tau_Aer, &
                                     wo_Aer, &
                                     g_Aer, &
                                     Tau_Ray, &
                                     Tau_gas, &
                                     scatangle, &
                                     Coszen, &
                                     Cossolzen, &
                                     Cloud_Albedo_View, &
                                     Cloud_Albedo_Sun, &
                                     Ref_ss)

   real, intent(in):: Tau_Aer
   real, intent(in):: wo_Aer
   real, intent(in):: g_Aer
   real, intent(in):: Tau_Ray
   real, intent(in):: Tau_Gas
   real, intent(in):: Scatangle
   real, intent(in):: Coszen
   real, intent(in):: Cossolzen
   real, intent(in):: Cloud_Albedo_View
   real, intent(in):: Cloud_Albedo_Sun
   real, intent(out):: Ref_ss

   real:: Airmass
   real:: P_Aer
   real:: P_Ray
   real:: Tau_Total
   real:: Tau_Scat_Total
   real:: Trans_Total
   real:: Tau_Iso_Total
   real:: Trans_Iso_Total_View
   real:: Trans_Iso_Total_Sun
   real:: Tau_Iso_Scat_Total
   real:: mu
   real:: Pf
   real:: wo
   real:: Ref_ss_a
   real:: Ref_ss_b
   real:: Ref_ss_c

      !--- compute cosine of scattering angle
      mu = cos(scatangle*dtor)

      !-- compute Rayleigh phase function
      Airmass = 1.0/Coszen + 1.0/Cossolzen
      P_Ray = 0.75*(1.0 + mu**2)

      !--- compute total transmission
      Tau_Total = Tau_Aer + Tau_Ray + Tau_gas
      Trans_Total = exp(-Tau_Total*Airmass)

      Tau_Iso_Total = (1.0-g_Aer)*Tau_Aer + Tau_Ray + Tau_gas
      Trans_Iso_Total_View = exp(-Tau_Iso_Total/Coszen)
      Trans_Iso_Total_Sun = exp(-Tau_Iso_Total/Cossolzen)

      !--- compute total scattering optical depth
      Tau_Scat_Total = wo_Aer*Tau_Aer + Tau_Ray
      Tau_Iso_Scat_Total = wo_Aer*(1.0-g_Aer)*Tau_Aer + Tau_Ray

      !--- single scatter albedo
      wo = (wo_Aer*Tau_Aer + Tau_Ray)/ ( Tau_Total )

      !aerosol phase function (Henyey-Greenstein)
      P_Aer = (1.0 - g_Aer**2)/( (1.0 + g_Aer**2 - 2.0*g_Aer*mu)**(1.5) )

      !--- compute effective phase function
      Pf = P_aer
      if (Tau_Scat_Total > 0.0) then
        Pf = (wo_Aer*Tau_Aer*P_Aer + Tau_Ray*P_Ray)/(Tau_Scat_Total)
      endif

      !--- compute single scatter reflectance (0-100%)

      !--- single scattering in the layer without hitting surface
      Ref_ss_a = wo*Pf/(4.0*Airmass*Coszen*Cossolzen) * (1.0 - Trans_Total )

      !--- single scattering in the layer without hitting surface
      Ref_ss_b = (Tau_Iso_Scat_Total / (2.0*Cossolzen)) * &
                  Trans_Iso_Total_View * Cloud_Albedo_View

      Ref_ss_c = (Tau_Iso_Scat_Total / (2.0*Coszen)) * &
                  Trans_Iso_Total_Sun * Cloud_Albedo_Sun

      Ref_ss = 100.0*(Ref_ss_a + Ref_ss_b + Ref_ss_c)

end subroutine COMPUTE_CLEAR_SKY_SCATTER

end module cx_atmospheric_correction_vis_mod
