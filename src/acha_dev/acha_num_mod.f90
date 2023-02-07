!$Id: acha_module.f90 3876 2020-06-18 13:34:40Z yli $
module ACHA_NUM_MOD
!---------------------------------------------------------------------
!
!----------------------------------------------------------------------

  use ACHA_SERVICES_MOD, only : &
           real4, int1, int4, real8, dtor, &
           PLANCK_RAD_FAST, PLANCK_TEMP_FAST, &
           ACHA_SYMBOL_STRUCT, &
           ACHA_RTM_NWP_STRUCT, &
           INVERT_MATRIX, LOCATE

  use KDTREE2_MODULE

  implicit none

  public:: GENERIC_PROFILE_INTERPOLATION
  public:: KNOWING_P_COMPUTE_T_Z
  public:: KNOWING_T_COMPUTE_P_Z
  public:: KNOWING_T_COMPUTE_P_Z_BOTTOM_UP
  public:: KNOWING_Z_COMPUTE_T_P
  public:: OPTIMAL_ESTIMATION
  public:: MEAN_SMOOTH2
  public:: KD_TREE_INTERP_2pred
  public:: KD_TREE_INTERP_3pred
  public:: COMPUTE_MEDIAN_SEGMENT
  public:: FINDGEN
  private:: COMPUTE_MEDIAN
  public:: COMPUTE_TIME_HOURS_ACHA

  !--- include the non-system specific variables
  include 'include/acha_parameters.inc'

  type(ACHA_SYMBOL_STRUCT), private :: Symbol

  integer, public, dimension(:), allocatable, save:: Chan_Idx_y      !deallocate this
  integer, public:: Inver_Top_Level_RTM
  integer, public:: Inver_Base_Level_RTM
  real, public:: Inver_Top_Height
  real, public:: Inver_Base_Height
  real, public:: Inver_Strength

  real, private, PARAMETER:: MISSING_VALUE_REAL4 = -999.0
  integer(kind=int1), private, PARAMETER:: MISSING_VALUE_integer1 = -128_int1
  integer(kind=int4), private, PARAMETER:: MISSING_VALUE_integer4 = -999

  real, private, PARAMETER:: REAL_INFINITY = huge(MISSING_VALUE_REAL4)

  contains 

  !---------------------------------------------------------------------------
  ! print vector of the same length
  !----------------------------------------------------------------------------
  subroutine showvector(A,title)
      real (kind=real4), dimension(:), intent(in):: A
      character(len=*), intent(in):: title
      integer:: ii
      write (unit=6,fmt=*)
      write (unit=6,fmt=*) title
      do ii = 1,size(A)
        write (unit=6,fmt="(f10.7)") A(ii)
      enddo
   end subroutine showvector
  !---------------------------------------------------------------------------
  ! print two vectors of the same length
  !----------------------------------------------------------------------------
   subroutine showvectors(A,B,title)
      real (kind=real4), dimension(:), intent(in):: A,B
      character(len=*), intent(in):: title
      integer:: ii
      write (unit=6,fmt=*)
      write (unit=6,fmt=*) title
      do ii = 1,size(A)
        write (unit=6,fmt="(f10.7,f10.7)") A(ii), B(ii)
      enddo
   end subroutine showvectors
  !-----------------------------------------------------------------------
  ! print a matrix
  !-----------------------------------------------------------------------
   subroutine showmatrix(A,title)
      real (kind=real4),dimension(:,:), intent(in):: A
      character(len=*), intent(in):: title
      integer:: ii
      write (unit=6,fmt=*)
      write (unit=6,fmt=*) title
      do ii = 1,size(A,1)
        write (unit=6,fmt="(16es12.4)") A(ii,:)
      enddo
   end subroutine showmatrix
!-----------------------------------------------------------------
! InterpoLate within profiles knowing Z to determine above cloud
! radiative terms used in forward model
!-----------------------------------------------------------------
function GENERIC_PROFILE_INTERPOLATION(X_value,X_Profile,Y_Profile)  &
            result(Y_value)

     real, intent(in):: X_value
     real, dimension(:), intent(in):: X_Profile
     real, dimension(:), intent(in):: Y_Profile
     real:: Y_value

     integer:: Lev_Idx
     real:: dx
     integer:: nlevels

     nlevels = size(X_Profile)

     !--- interpoLate pressure profile
     call LOCATE(X_Profile,nlevels,X_value,Lev_Idx)
     Lev_Idx = max(1,min(nlevels-1,Lev_Idx))

     dx = X_Profile(Lev_Idx+1) - X_Profile(Lev_Idx)

     !--- perform interpoLation
     if (dx /= 0.0) then
        Y_value = Y_Profile(Lev_Idx) +  &
                 (X_value - X_Profile(Lev_Idx))  * &
                 (Y_Profile(Lev_Idx+1) - Y_Profile(Lev_Idx)) / dx
     else
          Y_value = Y_Profile(Lev_Idx)
     endif

end function GENERIC_PROFILE_INTERPOLATION

!-----------------------------------------------------------------
! InterpoLate within profiles knowing P to determine T and Z
!-----------------------------------------------------------------
subroutine KNOWING_P_COMPUTE_T_Z(ACHA_RTM_NWP,P,T,Z,Lev_Idx)

     type(acha_rtm_nwp_struct), intent(in) :: ACHA_RTM_NWP
     real, intent(in):: P
     real, intent(out):: T
     real, intent(out):: Z
     integer, intent(out):: Lev_Idx
     real:: dp
     real:: dt
     real:: dz

     !--- initialize
     T = MISSING_VALUE_REAL4
     Z = MISSING_VALUE_REAL4
     Lev_Idx = MISSING_VALUE_integer4

     !--- check for missing
     if (P == MISSING_VALUE_REAL4) return

     !--- interpoLate pressure profile
     call LOCATE(ACHA_RTM_NWP%P_Prof,Num_Levels_RTM_Prof,P,Lev_Idx)
     Lev_Idx = max(1,min(Num_Levels_RTM_Prof-1,Lev_Idx))

     dp = ACHA_RTM_NWP%P_Prof(Lev_Idx+1) - ACHA_RTM_NWP%P_Prof(Lev_Idx)
     dt = ACHA_RTM_NWP%T_Prof(Lev_Idx+1) - ACHA_RTM_NWP%T_Prof(Lev_Idx)
     dz = ACHA_RTM_NWP%Z_Prof(Lev_Idx+1) - ACHA_RTM_NWP%Z_Prof(Lev_Idx)

     !--- perform interpoLation
       if (dp /= 0.0) then
           T = ACHA_RTM_NWP%T_Prof(Lev_Idx) + dt/dp * (P - ACHA_RTM_NWP%P_Prof(Lev_Idx))
           Z = ACHA_RTM_NWP%Z_Prof(Lev_Idx) + dz/dp * (P - ACHA_RTM_NWP%P_Prof(Lev_Idx))
       else
           T = ACHA_RTM_NWP%T_Prof(Lev_Idx)
           Z = ACHA_RTM_NWP%Z_Prof(Lev_Idx)
       endif

       !--- Some negative cloud heights are observed because  of bad height
       !--- NWP profiles.
       if (Z < 0) then
         Z = ZC_FLOOR
       endif

end subroutine KNOWING_P_COMPUTE_T_Z

!-----------------------------------------------------------------
! InterpoLate within profiles knowing Z to determine T and P
!-----------------------------------------------------------------
subroutine KNOWING_Z_COMPUTE_T_P(ACHA_RTM_NWP,P,T,Z,Lev_Idx)

     type(acha_rtm_nwp_struct), intent(in) :: ACHA_RTM_NWP
     real, intent(in):: Z
     real, intent(out):: T
     real, intent(out):: P
     integer, intent(out):: Lev_Idx
     real:: dp
     real:: dt
     real:: dz

     !--- initialize
     T = MISSING_VALUE_REAL4
     P = MISSING_VALUE_REAL4
     Lev_Idx = MISSING_VALUE_integer4

     !--- check for missing
     if (Z == MISSING_VALUE_REAL4) return

     !--- interpoLate pressure profile
     call LOCATE(ACHA_RTM_NWP%Z_Prof,Num_Levels_RTM_Prof,Z,Lev_Idx)
     Lev_Idx = max(1,min(Num_Levels_RTM_Prof-1,Lev_Idx))

     dp = ACHA_RTM_NWP%P_Prof(Lev_Idx+1) - ACHA_RTM_NWP%P_Prof(Lev_Idx)
     dt = ACHA_RTM_NWP%T_Prof(Lev_Idx+1) - ACHA_RTM_NWP%T_Prof(Lev_Idx)
     dz = ACHA_RTM_NWP%Z_Prof(Lev_Idx+1) - ACHA_RTM_NWP%Z_Prof(Lev_Idx)

     !--- perform interpoLation
     if (dz /= 0.0) then
           T = ACHA_RTM_NWP%T_Prof(Lev_Idx) + dt/dz * (Z - ACHA_RTM_NWP%Z_Prof(Lev_Idx))
           P = ACHA_RTM_NWP%P_Prof(Lev_Idx) + dp/dz * (Z - ACHA_RTM_NWP%Z_Prof(Lev_Idx))
     else
           T = ACHA_RTM_NWP%T_Prof(Lev_Idx)
           P = ACHA_RTM_NWP%P_Prof(Lev_Idx)
     endif

end subroutine KNOWING_Z_COMPUTE_T_P

!-----------------------------------------------------------------
! InterpoLate within profiles knowing T to determine P and Z
!-----------------------------------------------------------------
subroutine KNOWING_T_COMPUTE_P_Z(ACHA_RTM_NWP,Symbol, Cloud_Type, &
           P,T,Z,T_Tropo,Z_Tropo,P_Tropo,klev,ierr,Level_Within_Inversion_Flag)

     type(acha_rtm_nwp_struct), intent(in) :: ACHA_RTM_NWP
     type(acha_symbol_struct), intent(in) :: Symbol
     integer (kind=int1), intent(in):: Cloud_Type
     real, intent(in):: T
     real, intent(out):: P
     real, intent(out):: Z
     real, intent(in):: T_Tropo
     real, intent(in):: Z_Tropo
     real, intent(in):: P_Tropo
     integer, intent(out):: klev
     integer, intent(out):: ierr
     real:: dp
     real:: dt
     real:: dz
     integer:: kstart
     integer:: kend
     integer:: nlevels_temp
     integer, intent(out):: Level_Within_Inversion_Flag
    

     !--- initialization
     ierr = 0
     Z = MISSING_VALUE_REAL4
     P = MISSING_VALUE_REAL4
     klev = MISSING_VALUE_integer4

     !--- check for missing
     if (T == MISSING_VALUE_REAL4) return

     !--- test for existence of a valid solution with troposphere
     kstart = ACHA_RTM_NWP%Tropo_Level
     kend = ACHA_RTM_NWP%Sfc_Level
     Nlevels_Temp = kend - kstart + 1

     !--- check to see if warmer than max, than assume at surface
     if (T > maxval(ACHA_RTM_NWP%T_Prof(kstart:kend))) then
         P = ACHA_RTM_NWP%P_Prof(kend)
         Z = ACHA_RTM_NWP%Z_Prof(kend)
         klev = kend - 1
         ierr = 0
         return
     endif

     !--- check to see if colder than min, than assume above tropopause
     !--- and either limit height to tropopause or extrapoLate in stratosphere
     if (T < minval(ACHA_RTM_NWP%T_Prof(kstart:kend)) .or. T < T_Tropo) then
         if (ALLOW_STRATOSPHERE_SOLUTION_FLAG == 1 .and. Cloud_Type == Symbol%OVERSHOOTING_TYPE) then
           Z = Z_Tropo + (T - T_Tropo) / Dt_Dz_Strato
           P = P_Tropo + (Z - Z_Tropo) * Dp_Dz_Strato
         else
           P = P_Tropo
           Z = Z_Tropo
           klev = kstart + 1
         endif
         ierr = 0
         return
     endif

     !--- if there is an inversion, look below first
     Level_Within_Inversion_Flag = 0
     if (Inver_Top_Level_RTM > 0 .and. Inver_Base_Level_RTM > 0) then
         kstart = Inver_Top_Level_RTM
         kend =  Inver_Base_Level_RTM
         nlevels_temp = kend - kstart + 1
         call LOCATE(ACHA_RTM_NWP%T_Prof(kstart:kend),nlevels_temp,T,klev)
         if ((klev > 0) .and. (klev < nlevels_temp -1)) then
              klev = klev + kstart - 1
              Level_Within_Inversion_Flag = 1
         endif
      endif

    !--- if no solution within an inversion
    if (Level_Within_Inversion_Flag == 0) then
        kstart = ACHA_RTM_NWP%Tropo_Level
        kend = ACHA_RTM_NWP%Sfc_Level
        nlevels_temp = kend - kstart + 1
        call LOCATE(ACHA_RTM_NWP%T_Prof(kstart:kend),nlevels_temp,T,klev)
        if (klev == 0 .or. klev == nlevels_temp) then
        !if (klev == 0) then
            klev = minloc(abs(T-ACHA_RTM_NWP%T_Prof(kstart:kend)),1)
        endif
        klev = klev + kstart - 1
        klev = max(1,min(Num_Levels_RTM_Prof-1,klev))
    endif

    !--- General Inversion
    dp = ACHA_RTM_NWP%P_Prof(klev+1) - ACHA_RTM_NWP%P_Prof(klev)
    dt = ACHA_RTM_NWP%T_Prof(klev+1) - ACHA_RTM_NWP%T_Prof(klev)
    dz = ACHA_RTM_NWP%Z_Prof(klev+1) - ACHA_RTM_NWP%Z_Prof(klev)

    if (dt /= 0.0) then
        P = ACHA_RTM_NWP%P_Prof(klev) + dp/dt*(T-ACHA_RTM_NWP%T_Prof(klev))
        Z = ACHA_RTM_NWP%Z_Prof(klev) + dz/dt*(T-ACHA_RTM_NWP%T_Prof(klev))
    else
        P = ACHA_RTM_NWP%P_Prof(klev) 
        Z = ACHA_RTM_NWP%Z_Prof(klev)
    endif

end subroutine KNOWING_T_COMPUTE_P_Z

!-----------------------------------------------------------------
! InterpoLate within profiles knowing T to determine P and Z
! look at the bottom first and move up
!-----------------------------------------------------------------
subroutine KNOWING_T_COMPUTE_P_Z_BOTTOM_UP(ACHA_RTM_NWP, Symbol, &
      Cloud_Type,P,T,Z,T_Tropo,Z_Tropo,P_Tropo,klev,ierr,Level_Within_Inversion_Flag)

     type(acha_rtm_nwp_struct), intent(in) :: ACHA_RTM_NWP
     type(acha_symbol_struct), intent(in) :: Symbol
     integer (kind=int1), intent(in):: Cloud_Type
     real, intent(in):: T
     real, intent(out):: P
     real, intent(out):: Z
     real, intent(in):: T_Tropo
     real, intent(in):: Z_Tropo
     real, intent(in):: P_Tropo
     integer, intent(out):: klev
     integer, intent(out):: ierr 
     integer, intent(out):: Level_Within_Inversion_Flag

     real:: dp
     real:: dt
     real:: dz
     integer:: kstart
     integer:: kend
     integer:: nlevels_temp
     integer:: interp

     !--- initialization
     ierr = 0
     Z = MISSING_VALUE_REAL4
     P = MISSING_VALUE_REAL4
     klev = MISSING_VALUE_integer4
     Level_Within_Inversion_Flag = 0

     !--- check for missing
     if (T == MISSING_VALUE_REAL4) return

     !--- test for existence of a valid solution with troposphere
     kstart = ACHA_RTM_NWP%Tropo_Level
     kend = ACHA_RTM_NWP%Sfc_Level
     Nlevels_Temp = kend - kstart + 1

     !--- check to see if warmer than max, than assume at surface
     if (T > maxval(ACHA_RTM_NWP%T_Prof(kstart:kend))) then
         P = ACHA_RTM_NWP%P_Prof(kend)
         Z = ACHA_RTM_NWP%Z_Prof(kend)
         klev = kend - 1
         ierr = 0
         return
     endif

     !--- check to see if colder than min, than assume above tropopause
     !--- and either limit height to tropopause or extrapoLate in stratosphere
     if (T < minval(ACHA_RTM_NWP%T_Prof(kstart:kend)) .or. T < T_Tropo) then
         if (ALLOW_STRATOSPHERE_SOLUTION_FLAG == 1 .and. Cloud_Type == Symbol%OVERSHOOTING_TYPE) then
           Z = Z_Tropo + (T - T_Tropo) / Dt_Dz_Strato
           P = P_Tropo + (Z - Z_Tropo) * Dp_Dz_Strato
         else
           P = P_Tropo
           Z = Z_Tropo
           klev = kstart + 1
         endif
         ierr = 0
!        print *, "strat solution"
         return
     endif

    !--- find solution by looking at layers
    ierr = 1
    Level_Within_Inversion_Flag = 0
    do klev = kend, kstart - 1, -1

     !  print *, "test ", klev, T, ACHA_RTM_NWP%T_Prof(klev),ACHA_RTM_NWP%T_Prof(klev-1)

       if (abs(T - ACHA_RTM_NWP%T_Prof(klev)) < 0.1) then
               ierr = 0
               interp = 0
               exit
       endif

       if ((T >= ACHA_RTM_NWP%T_Prof(klev) .and. T < ACHA_RTM_NWP%T_Prof(klev-1)) .or. &
           (T <= ACHA_RTM_NWP%T_Prof(klev) .and. T > ACHA_RTM_NWP%T_Prof(klev-1))) then

           ierr = 0
           interp = 1

           if (T >= ACHA_RTM_NWP%T_Prof(klev) .and. T < ACHA_RTM_NWP%T_Prof(klev-1)) then
              Level_Within_Inversion_Flag = 1
           endif

           exit

       endif

    enddo

    if (ierr /= 0) then 
       return
    endif

    !--- General Inversion
    P = ACHA_RTM_NWP%P_Prof(klev) 
    Z = ACHA_RTM_NWP%Z_Prof(klev)

    dt = ACHA_RTM_NWP%T_Prof(klev) - ACHA_RTM_NWP%T_Prof(klev-1)

    if (interp == 1 .and. dt /= 0.0) then
      dp = ACHA_RTM_NWP%P_Prof(klev) - ACHA_RTM_NWP%P_Prof(klev-1)
      dz = ACHA_RTM_NWP%Z_Prof(klev) - ACHA_RTM_NWP%Z_Prof(klev-1)
      P = ACHA_RTM_NWP%P_Prof(klev) + dp/dt*(T-ACHA_RTM_NWP%T_Prof(klev))
      Z = ACHA_RTM_NWP%Z_Prof(klev) + dz/dt*(T-ACHA_RTM_NWP%T_Prof(klev))
    endif

!   print *, 'TEST ', T,ACHA_RTM_NWP%T_Prof(klev) , ACHA_RTM_NWP%T_Prof(klev-1), P, ACHA_RTM_NWP%P_Prof(klev) , ACHA_RTM_NWP%P_Prof(klev-1), interp, ierr

!   print *, "final T ", T, ACHA_RTM_NWP%T_Prof(klev), ACHA_RTM_NWP%T_Prof(klev-1)
!   print *, "final P ", P, ACHA_RTM_NWP%P_Prof(klev), ACHA_RTM_NWP%P_Prof(klev-1)

end subroutine KNOWING_T_COMPUTE_P_Z_BOTTOM_UP

!------------------------------------------------------------------------
! subroutine to compute the Iteration in x due to optimal
! estimation
!
! The notation in this routine follows that of Clive Rodgers (1976,2000)
!
! input to this routine:
! Iter_Idx - the number of the current Iteration
! Iter_Idx_Max - the maximum number of Iterations allowed
! nx - the number of x values
! ny - the number of y values
! Convergence_Criteria - the convergence criteria
! y - the vector of observations
! f - the vector of observations predicted by the forward model
! x - the vector of retrieved parameters
! x_Ap - the vector of the apriori estimate of the retrieved parameters
! K - the Kernel Matrix
! Sy - the covariance matrix of y and f
! Sa_inv - the inverse of the covariance matrix of x_Ap
! Delta_X_Max - the maximum step allowed for each Delta_X value
!
! output of this routine:
! Sx - the covariance matrix of x 
! Delta_X - the increment in x for the next Iteration
! Converged_Flag - flag indicating if convergence was met (yes or no)
! Fail_Flag - flag indicating if this process failed (yes or no)
!
! local variables:
! Sx_inv - the inverse of Sx
! Delta_X_dir - the unit direction vectors for delta-x 
! Delta_X_distance - the total length in x-space of Delta_X
! Delta_X_constrained - the values of Delta_X after being constrained
!-----------------------------------------------------------------------
subroutine OPTIMAL_ESTIMATION(Iter_Idx,Iter_Idx_Min,Iter_Idx_Max,nx,ny, &
                              Convergence_Criteria,Delta_X_Max, &
                              y,f,x,x_Ap,K,Sy,Sa_inv, &
                              Sx,AKM,Delta_x,Delta_x_prev, &
                              Conv_Test,Cost,Goodness,Converged_Flag,Fail_Flag)

  integer, intent(in):: Iter_Idx
  integer, intent(in):: Iter_Idx_Min
  integer, intent(in):: Iter_Idx_Max
  integer, intent(in):: ny
  integer, intent(in):: nx
  real(kind=real4), intent(in):: Convergence_Criteria
  real(kind=real4), dimension(:), intent(in):: Delta_X_Max
  real(kind=real4), dimension(:), intent(in):: y
  real(kind=real4), dimension(:), intent(in):: f
  real(kind=real4), dimension(:), intent(in):: x
  real(kind=real4), dimension(:), intent(in):: x_Ap
  real(kind=real4), dimension(:,:), intent(in):: K
  real(kind=real4), dimension(:,:), intent(in):: Sy
  real(kind=real4), dimension(:,:), intent(in):: Sa_inv
  real(kind=real4), dimension(:), intent(in):: Delta_x_prev
  real(kind=real4), dimension(:,:), intent(out):: Sx
  real(kind=real4), dimension(:,:), intent(out):: AKM
  real(kind=real4), intent(out):: Conv_Test
  real(kind=real4), intent(out):: Cost
  real(kind=real4), intent(out):: Goodness
  real(kind=real4), dimension(:), intent(out):: Delta_x
  real(kind=real4), dimension(ny,ny):: Sy_inv
  real(kind=real4), dimension(nx,nx):: Sx_inv
  real(kind=real4), dimension(nx):: Delta_x_constrained
  integer, intent(out):: Fail_Flag
  integer, intent(out):: Converged_Flag
  integer:: Singular_Flag
  integer:: ix
  integer:: m
  integer:: p
  logical :: Singular_Warned_Before
  integer:: Nbad
  
  Singular_Warned_Before = .false.

  m = size(Sy,1)
  p = size(Sx,1)

  Converged_Flag = 0
  Fail_Flag = 0
  Delta_X = MISSING_VALUE_REAL4
  Sx = MISSING_VALUE_REAL4


  ! + count(Sy < huge(Sy(0,0))--- check Sy
  Nbad = count(Sy > REAL_INFINITY) + count(Sy < -1.0*REAL_INFINITY)
  if (Nbad > 0) then
   print *, "Cloud Height warning ==> Bad Sy in ACHA "
   Fail_Flag = 1
   Converged_Flag = 0
   return 
  endif

  Singular_Flag =  INVERT_MATRIX(Sy, Sy_Inv, m)
  if (Singular_Flag == 1) then
   print *, "Cloud Height warning ==> Singular Sy in ACHA "
   Fail_Flag = 1
   Converged_Flag = 0
   return 
  endif

  !---- compute next step
  AKM = matmul(Transpose(K),matmul(Sy_inv,K))   !step saving
  Sx_inv = Sa_inv + AKM !(Eq.102 Rodgers)
  Singular_Flag =  INVERT_MATRIX(Sx_inv, Sx, p)
  if (Singular_Flag == 1 .and. .not.  Singular_Warned_Before ) then
   print *, "Cloud Height warning ==> Singular Sx in ACHA "
   print *, "x = ", x
   print *, "xa = ", x_Ap
   print *, "y = ", y
   print *, "f = ", f
   print *, "Sx_inv = ", Sx_Inv
   print *, "Sa_Inv = ", Sa_Inv
   print *, "AKM = ", AKM
   print *, "Sy = ", Sy
   print *, "Sy_Inv = ", Sy_Inv
   print *, "K  = ", K
   Converged_Flag = 0
   Fail_Flag = 1
   Singular_Warned_Before = .true.
   stop
   return
  endif

  Delta_x = matmul(Sx,(matmul(Transpose(K),matmul(Sy_inv,(y-f))) +  &
                       matmul(Sa_inv,x_Ap-x) ))

  !--------------------------------------------------------------
  ! compute averaging kernel matrix (note partialy computed above)
  !--------------------------------------------------------------
  AKM = matmul(Sx,AKM) 

  !--------------------------------------------------------------
  ! check for convergence
  !--------------------------------------------------------------

  !--- compute convergence metric
  Conv_Test = abs(sum(Delta_X*matmul(Sx_inv,Delta_X)))
  Goodness = sum((y-f)*matmul(Sy_Inv,y-f))
  Cost = sum((x-x_ap)*matmul(Sa_Inv,x-x_ap)) + Goodness

  !-------------------------------------------------------------------
  ! a direct constraint to avoid too large steps
  !-------------------------------------------------------------------
  do ix = 1,nx
     Delta_X_Constrained(ix) = sign( min( abs(Delta_X(ix)), Delta_X_Max(ix) ), Delta_X(ix) ) 
  enddo
  Delta_X  = Delta_X_Constrained

  ! if current and previous iteration Delta_x has opposite signs, reduce current
  ! magnitude
  if (Delta_x_prev(1) /= MISSING_VALUE_REAL4 .and. (Delta_x_prev(1)*Delta_x(1) <0 )) then
     Delta_X  = Delta_X_Constrained/5.0
  endif

  !--- check for exceeding allowed number of interactions
  if (Iter_Idx > Iter_Idx_Max) then
      Converged_Flag = 0
      Fail_Flag = 1
  endif

  !--- check for traditional convergence
  if (Conv_Test < Convergence_Criteria .and. Iter_Idx > Iter_Idx_Min) then
      Converged_Flag = 1
      Fail_Flag = 0
  endif

  if (Fail_Flag == 1 .and. Conv_Test < Convergence_Criteria) then
          print *, "HOLD ON ", Fail_Flag, Conv_Test, Convergence_Criteria, Iter_Idx
          stop
  endif

 end subroutine OPTIMAL_ESTIMATION
  !-------------------------------------------------------------------------------------------------
  ! Smooth a field using a mean over an area
  !
  ! Description
  !    Values of Z_in with Mask_In = 1 are used to populate pixels with Mask_Out
  !    = 1
  !    Z_Out is computed as the mean of Z_In*Mask_In over a box whose size is
  !    defined by N.
  !
  ! Input
  !    Mask_In - binary mask of point used as the source of the smoothing
  !    Mask_Out - binary mask of points to have a result upon exit
  !    Missin = missing value used as fill for Z_Out
  !    Count_Thresh - number of source points to compue an output
  !    N - half-width of smoothing box (x and y)
  !    Num_Elements = size of array in x-direction
  !    Num_Lines = size of array in y-direction
  !    Z_In - source values
  !    Z_Out - output values
  !    di = number of pixels to skip in the i direction (0=none,1=every other
  !    ...)
  !    dj = number of pixels to skip in the j direction (0=none,1=every other
  !    ...)
  !
  !-------------------------------------------------------------------------------------------------
  subroutine MEAN_SMOOTH2(Mask_In,Mask_Out,Missing,di,dj,N,Num_Elements, Num_Lines, Z_In,Z_Out)

     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_In
     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_Out
     real (kind=real4), intent(in):: Missing
     integer (kind=int4), intent(in):: di
     integer (kind=int4), intent(in):: dj

     integer (kind=int4), intent(in):: N
     integer (kind=int4), intent(in):: Num_Elements
     integer (kind=int4), intent(in):: Num_lines
     real (kind=real4), intent(in), dimension(:,:), target:: Z_In
     real (kind=real4), intent(out), dimension(:,:):: Z_Out
     integer (kind=int4), dimension(size(Mask_In,1),size(Mask_In,2)):: Count_Out
     integer:: i
     integer:: j
     real:: Count_Temporary
     real, pointer, dimension(:,:):: Z_In_Sub
     integer (kind=int1), pointer, dimension(:,:):: Mask_In_Sub
     integer (kind=int1), pointer, dimension(:,:):: Mask_Out_Sub
     integer:: i1,i2,j1,j2


     Z_Out = 0.0
     Count_Out = 0
     do j = 1 + dj, Num_Lines-dj, dj + 1

        j1 = min(Num_Lines,max(1,j - N))
        j2 = min(Num_Lines,max(1,j + N))

        do i = 1 + di, Num_Elements - di, di + 1

           if (Mask_In(i,j) == 0) cycle
           if (Z_out(i,j) > 0) cycle

           i1 = min(Num_Elements,max(1,i - N))
           i2 = min(Num_Elements,max(1,i + N))

           Mask_In_Sub => Mask_In(i1:i2,j1:j2)
           Mask_Out_Sub => Mask_Out(i1:i2,j1:j2)

           if (sum(Mask_Out_Sub) == 0) cycle

           Z_In_Sub => Z_In(i1:i2,j1:j2)
           Count_Temporary = sum(real(Mask_In_Sub))

           Z_Out(i1:i2,j1:j2) = Z_Out(i1:i2,j1:j2) + sum(Z_In_Sub*Mask_In_Sub) / Count_Temporary
           Count_Out(i1:i2,j1:j2) = Count_Out(i1:i2,j1:j2) + 1

           Mask_In_Sub => null()
           Mask_Out_Sub => null()
           Z_In_Sub => null()

        enddo
     enddo

     !--- make mean value
     where(Count_Out > 0)
         Z_Out = Z_Out / Count_Out
     endwhere
     where(Count_Out == 0)
         Z_Out = Missing
     endwhere

     !--- only values are missing where output mask is 0
     where(Mask_Out == 0)
       Z_Out = Missing
     endwhere

  end subroutine MEAN_SMOOTH2
!-----------------------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------------------
subroutine KD_TREE_INTERP_3pred(Mask_In,Mask_Out,pred_var1,pred_var2,pred_var3,Num_Elements, Num_Lines, nnbrute, Z_In,Z_Out)

     ! This function performs kd tree search. Z_In is the original data, and
     ! Z_Out is the output where indices flagged by Mask_Output are reassinged
     ! using average of values flagged by Mask_In.
     ! Currently values where Mask_In is set are not changed. nnbrute is the
     ! number of found closet indices for each search

     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_In
     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_Out
     integer (kind=int1), intent(in), dimension(:,:), target:: pred_var1
     real (kind=real4), intent(in), dimension(:,:), target:: pred_var2
     real (kind=real4), intent(in), dimension(:,:), target:: pred_var3
     integer (kind=int4), intent(in):: Num_Elements
     integer (kind=int4), intent(in):: Num_lines
     integer (kind=int4), intent(in):: nnbrute
     real (kind=real4), intent(in), dimension(:,:), target:: Z_In
     real (kind=real4), intent(out), dimension(:,:):: Z_Out

     ! define local variables, each KD-tree predictor is 1-D array
     type(kdtree2), pointer :: tree
     integer(kind=int4):: n_training, n_query
     integer:: i,j
     integer(kind=int4), dimension(:), allocatable:: ind_training, ind_query
     real(kind=real4), dimension(:), allocatable:: predictor_1
     real(kind=real4), dimension(:), allocatable:: predictor_2
     real(kind=real4), dimension(:), allocatable:: predictor_3
     real(kind=real4), dimension(:), allocatable:: out_temp_1d
     real(kind=real4), dimension(:), allocatable:: Z_Out_1d
     real(kind=real4), dimension(:,:), allocatable:: my_array
     real(kind=real4), dimension(:), allocatable:: query_vec
     type(kdtree2_result), dimension(:), allocatable :: results1

     integer:: Line_Idx
     integer:: Elem_Idx

     if (.not. allocated(predictor_1)) allocate(predictor_1(Num_Elements*Num_Lines))
     if (.not. allocated(predictor_2)) allocate(predictor_2(Num_Elements*Num_Lines))
     if (.not. allocated(predictor_3)) allocate(predictor_3(Num_Elements*Num_Lines))
     if (.not. allocated(out_temp_1d)) allocate(out_temp_1d(Num_Elements*Num_Lines))
     if (.not. allocated(Z_Out_1d)) allocate(Z_Out_1d(Num_Elements*Num_Lines))

     Z_Out_1d = MISSING_VALUE_REAL4

     ! find indices of training and query variables
     ! values at training indices (Mask_In) can be computed again but original
     ! values are kept here
     n_training = COUNT(Mask_In == 1)
     n_query = COUNT(Mask_Out == 1) - Count(Mask_Out == 1 .and. Mask_In == 1)

     ! perform kd-tree only if both training and query indices are found
     if (n_training > 0 .and. n_query > 0) then
        ! convert 2d to 1d array

        predictor_1 = reshape(real(pred_var1),(/Num_Elements*Num_Lines/))

        predictor_2 = reshape(pred_var2,(/Num_Elements*Num_Lines/))

        predictor_3 = reshape(pred_var3,(/Num_Elements*Num_Lines/))

        out_temp_1d = reshape(Z_In,(/Num_Elements*Num_Lines/))

        if (.not. allocated(ind_training)) allocate(ind_training(n_training))
        if (.not. allocated(ind_query)) allocate(ind_query(n_query))
        ! search for training and query indices and stored
        i = 1
        j = 1
        do Line_Idx = 1, Num_Lines
        do Elem_Idx = 1, Num_Elements
!          if (Mask_In(Elem_Idx,Line_Idx) == 1 .and.
!          pred_var1(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
          if (Mask_In(Elem_Idx,Line_Idx) == 1 .and. pred_var1(Elem_Idx,Line_Idx) /= real(MISSING_VALUE_integer1) &
              .and. pred_var2(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var3(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              ind_training(i) = Elem_Idx + Num_Elements*(Line_Idx-1)
              i = i+1
          endif
          if (Mask_Out(Elem_Idx,Line_Idx) == 1 .and. Mask_In(Elem_Idx,Line_Idx) == 0 &
              .and. pred_var1(Elem_Idx,Line_Idx) /= real(Missing_Value_integer1) &
!              .and. pred_var1(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var2(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var3(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              ind_query(j) = Elem_Idx + Num_Elements*(Line_Idx-1)
              j = j+1
          endif
          enddo
        enddo

        if (minval(out_temp_1d(ind_training)) == Missing_Value_Real4) print *, 'kdtree training is not correct'

        if (.not. allocated(my_array)) allocate(my_array(3,n_training))
        if (.not. allocated(query_vec)) allocate(query_vec(3))

        ! assign training data
        my_array(1,:) = predictor_1(ind_training)
        my_array(2,:) = predictor_2(ind_training)
        my_array(3,:) = predictor_3(ind_training)

        ! create tree, set sort to true slows it down but the output indices are
        ! ordered from closet to farthest; set rearrange as true speeds searches
        ! but requires extra memory

        tree => kdtree2_create(my_array,sort=.true.,rearrange=.true.)  ! this is how you create a tree.

        if (.not. allocated(results1)) allocate(results1(nnbrute))

        ! set 1d output variable values at training indices
        do i = 1,n_training
           Z_Out_1d(ind_training(i)) = out_temp_1d(ind_training(i))
        enddo

        ! perform tree search for each query index
        do i = 1,n_query
           query_vec(1) = predictor_1(ind_query(i))
           query_vec(2) = predictor_2(ind_query(i))
           query_vec(3) = predictor_3(ind_query(i))

           ! results1 has both indices and distances 
           call kdtree2_n_nearest(tp=tree, qv=query_vec, nn = nnbrute, results=results1)

           ! average values for the all found indices
           Z_Out_1d(ind_query(i)) = sum(out_temp_1d(ind_training(results1%idx)))/size(results1%idx)
        enddo

        ! destroy tree and release memory
        call kdtree2_destroy(tree)

     endif

     ! Change 1d array back to 2d; if no kdtree search is performed, array is
     ! empty
     Z_Out = reshape(Z_Out_1d,(/Num_Elements,Num_Lines/))


     if (allocated(predictor_1)) deallocate(predictor_1)
     if (allocated(predictor_2)) deallocate(predictor_2)
     if (allocated(predictor_3)) deallocate(predictor_3)
     if (allocated(out_temp_1d)) deallocate(out_temp_1d)
     if (allocated(Z_Out_1d)) deallocate(Z_Out_1d)
     if (allocated(query_vec)) deallocate(query_vec)
     if (allocated(my_array)) deallocate(my_array)
     if (allocated(ind_training)) deallocate(ind_training)
     if (allocated(ind_query)) deallocate(ind_query)
     if (allocated(results1)) deallocate(results1)

  end subroutine KD_TREE_INTERP_3pred

!-----------------------------------------------------------------------------------------------------
! 
!-----------------------------------------------------------------------------------------------------
subroutine KD_TREE_INTERP_2pred(Mask_In,Mask_Out,pred_var1,pred_var2,Num_Elements, Num_Lines, nnbrute, Z_In,Z_Out)

     ! This function performs kd tree search. Z_In is the original data, and
     ! Z_Out is the output where indices flagged by Mask_Output are reassinged
     ! using average of values flagged by Mask_In.
     ! Currently values where Mask_In is set are not changed. nnbrute is the
     ! number of found closet indices for each search

     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_In
     integer (kind=int1), intent(in), dimension(:,:), target:: Mask_Out
     real (kind=real4), intent(in), dimension(:,:), target:: pred_var1
     real (kind=real4), intent(in), dimension(:,:), target:: pred_var2
     integer (kind=int4), intent(in):: Num_Elements
     integer (kind=int4), intent(in):: Num_lines
     integer (kind=int4), intent(in):: nnbrute
     real (kind=real4), intent(in), dimension(:,:), target:: Z_In
     real (kind=real4), intent(out), dimension(:,:):: Z_Out

     ! define local variables, each KD-tree predictor is 1-D array
     type(kdtree2), pointer :: tree
     integer(kind=int4):: n_training, n_query
     integer:: i,j
     integer(kind=int4), dimension(:), allocatable:: ind_training, ind_query
     real(kind=real4), dimension(:), allocatable:: predictor_1
     real(kind=real4), dimension(:), allocatable:: predictor_2
     real(kind=real4), dimension(:), allocatable:: out_temp_1d
     real(kind=real4), dimension(:), allocatable:: Z_Out_1d
     real(kind=real4), dimension(:,:), allocatable:: my_array
     real(kind=real4), dimension(:), allocatable:: query_vec
     type(kdtree2_result), dimension(:), allocatable :: results1

     integer:: Line_Idx
     integer:: Elem_Idx

     if (.not. allocated(predictor_1)) allocate(predictor_1(Num_Elements*Num_Lines))
     if (.not. allocated(predictor_2)) allocate(predictor_2(Num_Elements*Num_Lines))
     if (.not. allocated(out_temp_1d)) allocate(out_temp_1d(Num_Elements*Num_Lines))
     if (.not. allocated(Z_Out_1d)) allocate(Z_Out_1d(Num_Elements*Num_Lines))

     Z_Out_1d = MISSING_VALUE_REAL4

     ! find indices of training and query variables
     ! values at training indices (Mask_In) can be computed again but original
     ! values are kept here
     n_training = COUNT(Mask_In == 1)
     n_query = COUNT(Mask_Out == 1) - Count(Mask_Out == 1 .and. Mask_In == 1)

     ! perform kd-tree only if both training and query indices are found
     if (n_training > 0 .and. n_query > 0) then
        ! convert 2d to 1d array

        predictor_1 = reshape(pred_var1,(/Num_Elements*Num_Lines/))

        predictor_2 = reshape(pred_var2,(/Num_Elements*Num_Lines/))

        out_temp_1d = reshape(Z_In,(/Num_Elements*Num_Lines/))

        if (.not. allocated(ind_training)) allocate(ind_training(n_training))
        if (.not. allocated(ind_query)) allocate(ind_query(n_query))
        ! search for training and query indices and stored
        i = 1
        j = 1
        do Line_Idx = 1, Num_Lines
        do Elem_Idx = 1, Num_Elements
          if (Mask_In(Elem_Idx,Line_Idx) == 1  &
              .and. pred_var1(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var2(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              ind_training(i) = Elem_Idx + Num_Elements*(Line_Idx-1)
              i = i+1
          endif
          if (Mask_Out(Elem_Idx,Line_Idx) == 1 .and. Mask_In(Elem_Idx,Line_Idx) == 0 &
              .and. pred_var1(Elem_Idx,Line_Idx) /= Missing_Value_Real4 &
              .and. pred_var2(Elem_Idx,Line_Idx) /= Missing_Value_Real4) then
              ind_query(j) = Elem_Idx + Num_Elements*(Line_Idx-1)
              j = j+1
          endif
          enddo
        enddo

        if (minval(out_temp_1d(ind_training)) == Missing_Value_Real4) print *, 'kdtree training is not correct'

        if (.not. allocated(my_array)) allocate(my_array(2,n_training))
        if (.not. allocated(query_vec)) allocate(query_vec(2))

        ! assign training data
        my_array(1,:) = predictor_1(ind_training)
        my_array(2,:) = predictor_2(ind_training)

        ! create tree, set sort to true slows it down but the output indices are
        ! ordered from closet to farthest; set rearrange as true speeds searches
        ! but requires extra memory

        tree => kdtree2_create(my_array,sort=.true.,rearrange=.true.)  ! this is how you create a tree.

        if (.not. allocated(results1)) allocate(results1(nnbrute))

        ! set 1d output variable values at training indices
        do i = 1,n_training
           Z_Out_1d(ind_training(i)) = out_temp_1d(ind_training(i))
        enddo

        ! perform tree search for each query index
        do i = 1,n_query
           query_vec(1) = predictor_1(ind_query(i))
           query_vec(2) = predictor_2(ind_query(i))

           ! results1 has both indices and distances 
           call kdtree2_n_nearest(tp=tree, qv=query_vec, nn = nnbrute, results=results1)

           ! average values for the all found indices
           Z_Out_1d(ind_query(i)) = sum(out_temp_1d(ind_training(results1%idx)))/size(results1%idx)
        enddo

        ! destroy tree and release memory
        call kdtree2_destroy(tree)

     endif

     ! Change 1d array back to 2d; if no kdtree search is performed, array is
     ! empty
     Z_Out = reshape(Z_Out_1d,(/Num_Elements,Num_Lines/))


     if (allocated(predictor_1)) deallocate(predictor_1)
     if (allocated(predictor_2)) deallocate(predictor_2)
     if (allocated(out_temp_1d)) deallocate(out_temp_1d)
     if (allocated(Z_Out_1d)) deallocate(Z_Out_1d)
     if (allocated(query_vec)) deallocate(query_vec)
     if (allocated(my_array)) deallocate(my_array)
     if (allocated(ind_training)) deallocate(ind_training)
     if (allocated(ind_query)) deallocate(ind_query)
     if (allocated(results1)) deallocate(results1)

  end subroutine KD_TREE_INTERP_2pred

 !----------------------------------------------------------------------
 ! subroutine COMPUTE_MEDIAN_SEGMENT(z,mask,n,imin,imax,jmin,jmax,
 !                                   z_median,z_std_median)
 !
 ! Compute standard deviaion of an array wrt to the median
 !----------------------------------------------------------------------
 subroutine COMPUTE_MEDIAN_SEGMENT(z,mask,n,nx,ny,z_median, z_std_median)
  real(kind=real4), dimension(:,:), intent(in):: z
  integer(kind=int1), dimension(:,:), intent(in):: mask
  real(kind=real4), dimension(:,:), intent(out):: z_median
  real(kind=real4), dimension(:,:), intent(out), optional:: z_std_median
  integer, intent(in):: n
  integer, intent(in):: nx
  integer, intent(in):: ny
  integer:: i
  integer:: j
  integer:: i1
  integer:: i2
  integer:: j1
  integer:: j2
  real(kind=real4) :: z_mean, z_std

  do i = 1, nx 
    do j = 1, ny

     j1 = max(1,j-n)   !top index of local array
     j2 = min(ny,j+n)   !bottom index of local array
     i1 = max(1,i-n)   !left index of local array
     i2 = min(nx,i+n)   !right index of local array

     !--- compute median
     call COMPUTE_MEDIAN(z(i1:i2,j1:j2),mask(i1:i2,j1:j2),z_median(i,j), &
                         z_mean,z_std)
     if (present(z_std_median)) z_std_median(i,j) = z_std

     enddo
  enddo

 end subroutine COMPUTE_MEDIAN_SEGMENT
 !==============================================================
 ! subroutine COMPUTE_MEDIAN(z,mask,z_median,z_mean,z_std_median)
 !
 ! Median filter
 !==============================================================
 subroutine COMPUTE_MEDIAN(z,mask,z_median,z_mean,z_std_median)

! The purpose of this function is to find
! median (emed), minimum (emin) and maximum (emax)

 ! for the array elem with nelem elements.

 real, dimension(:,:), intent(in):: z
 real, intent(out):: z_median
 real, intent(out):: z_mean
 real, intent(out):: z_std_median
 integer(kind=int1), dimension(:,:), intent(in):: mask
 integer:: i,j,k,nx,ny,nelem
 real, dimension(:), allocatable::x
 real(kind=real4):: u

 z_median = missing_value_real4
 z_std_median = missing_value_real4
 z_mean = missing_value_real4

 nx = size(z,1)
 ny = size(z,2)

 nelem = nx * ny

 allocate(x(nelem))
 x = 0.0

 k = 0
 do i = 1, nx
   do j = 1, ny
      if (mask(i,j) == 0 .and. z(i,j) /= missing_value_real4) then
           k = k + 1
           x(k) = z(i,j)
      endif
  enddo
 enddo

 nelem = k

 if (nelem < 1) then
     if (allocated(x)) deallocate(x)
     return
 endif
!--- sort the array into ascending order
  do i=1,nelem-1
   do j=i+1,nelem
    if(x(j)<x(i))then
     u=x(j)
     x(j)=x(i)
     x(i)=u
    end if
   end do
  end do

 !---- pick the median
  if(mod(nelem,2)==1)then
   i=nelem/2+1
   z_median=x(i)
  else
   i=nelem/2
   z_median=(x(i)+x(i+1))/2
   end if

 !--- compute standard deviation wrt median
  z_mean = sum(x(1:nelem))/nelem
  z_std_median = sqrt(sum((x(1:nelem) - z_median)**2) / nelem)


  if (allocated(x)) deallocate(x)

 end subroutine COMPUTE_MEDIAN

 !----------------------------------------------------------------------
 function findgen( n ) result( r )
    integer, intent(in) :: n
    integer:: i
    real, dimension(n) :: r
    do i = 1,n
        r(i) = i
    enddo
 end function findgen

 !---------------------------------------------------------------
 function COMPUTE_TIME_HOURS_ACHA()   result(time_hours)
   character(len=8):: system_date
   character(len=10):: system_time
   character(len=5):: system_time_zone
   integer, dimension(8):: system_time_value

   real:: time_hours

   call DATE_AND_TIME(system_date,system_time,system_time_zone, system_time_value)

   time_hours = real(system_time_value(5)) +  &
                     (real(system_time_value(6)) + &
                      real(system_time_value(7) + &
                      real(system_time_value(8))/1000.0)/60.0)/60.0
   return

  end function COMPUTE_TIME_HOURS_ACHA

!----------------------------------------------------------------------
! End of Module
!----------------------------------------------------------------------

end module ACHA_NUM_MOD
