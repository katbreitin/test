! $Id: cx_science_tools_mod.f90 4061 2021-01-09 01:31:15Z heidinger $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE
!
! NAME: cx_science_tools_mod.f90 (src)
!       CX_SCIENCE_TOOLS_MOD (program)
!
! PURPOSE: library of useful scientifix routines and  functions
!
! Description: 
!
! AUTHORS:
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
! public:: 
!   LOCATE
!   JULIAN
!   COMPUTE_MONTH
!   COMPUTE_DAY
!   VAPOR
!   VAPOR_ICE
!   NOIST_ADIABATIC_LAPSE_RATE
!
!--------------------------------------------------------------------------------------
 module CX_SCIENCE_TOOLS_MOD
  use CONSTANTS_MOD
  
  implicit none
  private
  public::  &
           VAPOR, &
           VAPOR_ICE, &
           MOIST_ADIABATIC_LAPSE_RATE, &
           WIND_SPEED, &
           WIND_DIRECTION, &
           COMPUTE_LFC_EL_HEIGHT
  contains



!----------------------------------------------------------------
! compute level of free convection and equilibrium level
!
! input:
! P = Pressure Profile (hPa)
! T = Temperature Profile (K)
! Td_Sfc = Dew Point Temperature at Surface (K)
! Z = height profile (m)
!
! output:
! LFC = level of free convection (m)
! EL = equilibrium level (m)
!----------------------------------------------------------------
subroutine COMPUTE_LFC_EL_HEIGHT(Sfc_Level,Tropo_Level,Nlevels, &
                                 P,T,Td_Sfc,Z,LFC,EL)
   integer, intent(in):: Nlevels
   integer(kind=int1),intent(in):: Sfc_Level, Tropo_Level
   real, intent(in):: Td_Sfc
   real, dimension(:), intent(in):: P,T,Z
   real, intent(out):: LFC, EL
   real:: malr
   real, dimension(nlevels):: T_Moist
   integer, dimension(nlevels):: xMask
   integer:: Lev_Idx
   integer, dimension(1):: Idx

   T_Moist = MISSING_VALUE_REAL4
   EL = MISSING_VALUE_REAL4
   LFC = MISSING_VALUE_REAL4

   T_Moist(Sfc_Level:Nlevels) = Td_Sfc

   do Lev_Idx = Sfc_Level,Tropo_Level, -1
     Malr = MOIST_ADIABATIC_LAPSE_RATE(T(Lev_Idx),P(Lev_Idx))  !K/km
     T_Moist(Lev_Idx-1) = T_Moist(Lev_Idx) -  &
                          Malr*(Z(Lev_Idx-1)-Z(Lev_Idx)) / 1000.0  !Z is m

     !--- level of free convection
     if ((LFC == MISSING_VALUE_REAL4) .and. &
         (T_Moist(Lev_Idx-1) > T(Lev_Idx-1)) .and. &
         (T_Moist(Lev_Idx) < T(Lev_Idx))) then

         LFC = Z(Lev_Idx-1)

     endif

     !--- find equilibrium level
     if ((EL == MISSING_VALUE_REAL4) .and. &
         (T_Moist(Lev_Idx-1) < T(Lev_Idx-1)) .and. &
         (T_Moist(Lev_Idx) > T(Lev_Idx))) then

         EL = Z(Lev_Idx-1)

     endif

   enddo

!---- findloc is in gfortran 9
!  xMask = 0
!  where (T_Moist > T)
!    xMask = 1
!  endwhere

!  Idx = findloc(xMask,1,back=.true.)
!  LFC = Z(Idx(1))

!  Idx = findloc(xMask,1)
!  EL = Z(Idx(1))
   
end subroutine COMPUTE_LFC_EL_HEIGHT

!----------------------------------------------------------------
! moist adiabatic lapse rate
!
! T = Temperature (K)
! P = Pressure (hPa)
!
! output
! malr = moist adiabatic lapse rate (K/km)
!
! method: http://www.theweatherprediction.com/habyhints/161/
!
!alternateL http://hogback.atmos.colostate.edu/group/dave/pdf/
!           Moist_adiabatic_lapse_rate.pdf
!----------------------------------------------------------------
 function MOIST_ADIABATIC_LAPSE_RATE(T,P) result(malr)

  real, intent(in):: T, P
  real:: es, ws, dws_dt, qp
  real:: Tp, esp, wsp
  real:: malr
  real, parameter:: eps = 0.622
  real, parameter:: dalr = 9.8
  real, parameter:: L = 2.453e06 !J/kg Latent heat of vaporization 
  real, parameter:: cp = 1004.0 ! J/kg-K

  if (T > 253.0) then
    es = VAPOR(T)
  else
    es = VAPOR_ICE(T)
  endif

  ws = 0.622*es / (P -es)

  Tp = T + 1

  if (Tp > 253.0) then
    esp = VAPOR(Tp)
  else
    esp = VAPOR_ICE(Tp)
  endif

  wsp = 0.622*esp / (P -esp)

  dws_dt = wsp - ws

  malr = dalr / (1.0 + L/Cp*dWs_dT)  

 return

 end function MOIST_ADIABATIC_LAPSE_RATE

!----------------------------------------------------------------
! functions to compute some needed water vapor parameters
!----------------------------------------------------------------
 function VAPOR(T) result(es)
                                                                     
!  T in Kelvin                                                          
!  es in mbar

  implicit none
  real, intent (in) :: T
  real :: es

   es = 6.112 * exp(17.67 * (T-273.16) / (T - 29.66))

  return 
end function VAPOR

!---- saturation vapor pressure for ice
function VAPOR_ICE(T) result(es)
   implicit none
   real, intent(in):: T
   real:: es
     es = 6.1078 * exp(21.8745584 * (T-273.16) / (T - 7.66))
  return 
end function VAPOR_ICE

!---------------------------------------------------------------------
! SUBPROGRAM:  W3FC05        EARTH U,V WIND COMPONENTS TO DIR AND SPD
!   PRGMMR: CHASE            ORG: NMC421      DATE:88-10-26
!
! ABSTRACT: GIVEN THE TRUE (EARTH ORIENTED) WIND COMPONENTS
!   COMPUTE THE WIND DIRECTION AND SPEED.
!   INPUT WINDS AT THE POLE ARE ASSUMED TO FOLLOW THE WMO
!   CONVENTIONS, WITH THE OUTPUT DIRECTION COMPUTED IN ACCORDANCE
!   WITH WMO STANDARDS FOR REPORTING WINDS AT THE POLE.
!   (SEE OFFICE NOTE 241 FOR WMO DEFINITION.)
!
! PROGRAM HISTORY LOG:
!   81-12-30  STACKPOLE, JOHN
!   88-10-19  CHASE, P.   ALLOW OUTPUT VALUES TO OVERLAY INPUT
!   89-01-21  R.E.JONES   CONVERT TO MICROSOFT FORTRAN 4.10
!   90-06-11  R.E.JONES   CONVERT TO SUN FORTRAN 1.3
!   91-03-30  R.E.JONES   SiliconGraphics FORTRAN
!
! USAGE:    call W3FC05 (U, V, DIR, SPD)
!
!   INPUT ARGUMENT LIST:
!     U        - real*4 EARTH-ORIENTED U-COMPONENT
!     V        - real*4 EARTH-ORIENTED V-COMPONENT
!
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!     DIR      - real*4 WIND DIRECTION, DEGREES.  VALUES WILL
!                BE FROM 0 TO 360 INCLUSIVE.
!     SPD      - real*4 WIND SPEED IN SAME UNITS AS INPUT
!---------------------------------------------------------------------
subroutine WIND_SPEED_AND_DIRECTION(u,v,dir,spd)
  real, intent(in) :: u
  real, intent(in) :: v
  real, intent(out) :: spd
  real, intent(out) :: dir
 

  real, parameter:: SPDTST = 1.0e-10
  real, parameter:: RTOD = 57.2957795
  real, parameter:: dchalf = 180.0

  spd = Missing_Value_Real4
  dir = Missing_Value_Real4
  if (u == Missing_Value_Real4 .or. v == Missing_Value_Real4) then
     return
  endif 
  spd = sqrt(u * u + v * v)
  if (spd < SPDTST) THEN
        dir = 0.0
  else 
        dir = atan2(u,v) * RTOD + DCHALF
  endif

  return

end subroutine WIND_SPEED_AND_DIRECTION
!-------------------------------------------------------------------------------------
! elemental funcions wind_speed and wind_direction taken from W3FC03 from UCAR
!
! W3FC05 header documentation follows:
!
! ABSTRACT: GIVEN THE TRUE (EARTH ORIENTED) WIND COMPONENTS
!   COMPUTE THE WIND DIRECTION AND SPEED.
!   INPUT WINDS AT THE POLE ARE ASSUMED TO FOLLOW THE WMO
!   CONVENTIONS, WITH THE OUTPUT DIRECTION COMPUTED IN ACCORDANCE
!   WITH WMO STANDARDS FOR REPORTING WINDS AT THE POLE.
!   (SEE OFFICE NOTE 241 FOR WMO DEFINITION.)
!
! PROGRAM HISTORY LOG:
!   81-12-30  STACKPOLE, JOHN
!   88-10-19  CHASE, P.   ALLOW OUTPUT VALUES TO OVERLAY INPUT
!   89-01-21  R.E.JONES   CONVERT TO MICROSOFT FORTRAN 4.10
!   90-06-11  R.E.JONES   CONVERT TO SUN FORTRAN 1.3
!   91-03-30  R.E.JONES   SiliconGraphics FORTRAN
!
! USAGE:    CALL W3FC05 (U, V, DIR, SPD)
!
!   INPUT ARGUMENT LIST:
!     U        - real*4 EARTH-ORIENTED U-COMPONENT
!     V        - real*4 EARTH-ORIENTED V-COMPONENT
!
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!     DIR      - real*4 WIND DIRECTION, DEGREES.  VALUES WILL
!                BE FROM 0 TO 360 INCLUSIVE.
!     SPD      - real*4 WIND SPEED IN SAME UNITS AS INPUT
!-------------------------------------------------------------------------------------
real elemental function WIND_SPEED ( u ,v )

  real, intent(in)::  u
  real, intent(in):: v

  if (u == Missing_Value_Real4 .or. v == Missing_Value_Real4) then
     wind_speed = Missing_Value_Real4
  else
     wind_speed = sqrt(u * u + v * v)
  endif

end function wind_speed
!-------------------------------------------------------------------------------------
! taken from W3FC03 from UCAR
!-------------------------------------------------------------------------------------
real elemental function WIND_DIRECTION ( u ,v )

  real, intent(in)::  u
  real, intent(in):: v
  real, parameter:: rtod = 57.2957795
  real, parameter:: dchalf = 180.0
  real, parameter:: spdtst = 1.0e-10

  if (u == Missing_Value_Real4 .or. v == Missing_Value_Real4) then
     wind_direction = Missing_Value_Real4
  else
     if (abs(u) < spdtst .and. abs(v) < spdtst) then
      wind_direction =  0.0
     else
      wind_direction = atan2(u,v) * rtod + dchalf
     endif
  endif

end function wind_direction

!------------------------------------------------------------------------------------- 
end module CX_SCIENCE_TOOLS_MOD

