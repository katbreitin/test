! $Id: muri_retrieval_mod.f90 2300 2017-08-29 10:33:06Z awalther $

module muri_land_retrieval_mod
  use muri_definitions_mod, only: &
  muri_input_type &
  , muri_output_type !&
  !		, surface_refl_type

  use muri_land_lut_mod, only: &
  land_lut

  use lib_array, only:interp1d
  use aw_lib_array, only: interp3d

  use univ_fp_comparison_mod, only: operator(.EQfp.)

  implicit none

  public :: muri_land_algorithm



contains

  !
  !
  !
  subroutine muri_land_algorithm (inp, out)
    implicit none


    type( muri_input_type), intent(in) :: inp
    type( muri_output_type),intent(out) :: out

    integer :: n_opt ,i_size
    integer, parameter :: N_FMR=11
    integer, parameter :: N_BANDS = 3
    integer, parameter :: i_band2p3=3 ! the 3rd index is 2p3 2.3 um wavelengths
    real :: fmr
    real, allocatable:: refl_toa (:,:,:)


    integer :: i_cha,i_opt, i_fmr
    real, allocatable :: refl_reference (:)
    real :: aot_temp(N_FMR)
    real :: refl_corrsp(N_FMR,3)
    real :: aot_allbands(N_FMR,3)
    real :: n1 (N_FMR)
    real :: n2
    real :: err(N_FMR,3)
    real :: err_640 (N_FMR)
    real :: opt_0470(8),opt_0640(8) ! NOPT
    real :: opt_b1,opt_b3
    real :: wavb1,wavb3

    real :: aod_nleta_temp
    real :: err_nleta_temp
    real :: fmf_nleta_temp
    integer :: idx(1)

    real :: val
    real :: aod_allbands(3)
    integer :: CHANNEL_REFERENCE

    real :: MVI
    real :: yint644,slope644
    real :: yint466,slope466

    real :: P1,P2
    real, allocatable :: RHO_surf_2p3(:,:)
    real :: RHO_surf_644,RHO_surf_466
    real :: denb,numb

    real :: lut_toa_ref(8,3,2)

    integer :: band(3)

    ! surface elev correction
    real :: MHGHT
    integer, parameter :: NLTAU=8 ! number of optical thickness
    integer, parameter :: NLWAV=3 ! number of wavelength/band
    integer, parameter :: NLSIZE=2 ! number of optical thickness
    REAL :: INT_NL(NLTAU,NLWAV,NLSIZE)
    REAL :: Fd_NL(NLTAU,NLWAV,NLSIZE), T_NL(NLTAU,NLWAV,NLSIZE)
    REAL :: SBAR_NL(NLTAU,NLWAV,NLSIZE),OPTH_NL(NLTAU,NLWAV,NLSIZE)
    REAL :: REF_RAY_NL(NLWAV)
    real :: WAV_NL(NLWAV),EQWAV_NL(NLWAV)


    wavb1=470
    wavb3=640

    WAV_NL=(/0.47, 0.64, 2.3/)
    band=(/1,3,6/)

    MHGHT=inp % surf_elev
    MHGHT=MHGHT/1000


    ! call inp%info
    CHANNEL_REFERENCE = 1  ! channel reference over land is 0.47


    err_nleta_temp = huge(err_nleta_temp)
    err_640 = huge(err_640)
    aod_nleta_temp=-99.00



    !********************************************************
    ! surface reflectance parameter to estimate
    ! this part can be modified to subroutine or function later

    !MVIB6B4=(inp % rfl(4)-inp % rfl(6))/(inp % rfl(4)+inp % rfl(6))
    !MVI=(0.23*MVIB6B4**2)+(0.67*MVIB6B4)+0.008

    MVI=(inp % rfl(4)-inp % rfl(6))/(inp % rfl(4)+inp % rfl(6))


    yint466=0.005
    slope466=0.49

    ! yint644=0.0 to  add yint644=0.025

    if(MVI.gt.0.75) then
      yint644=0.025
      slope644=0.58

    else if(MVI.lt.0.25) then
      yint644=0.025
      slope644=0.48

    else
      yint644=0.025
      slope644=0.48+0.2*(MVI-0.25)

    end if
    !print*, 'slopes 644 should be between 0.48 to 0.58 :',slope644


    yint644=yint644-0.00025*inp%scat_angle+0.033684
    slope644=slope644+0.002*inp%scat_angle-0.27 +0.025 ! increase slope 0.025

    ! print*,'scatteringh angle',inp%scat_angle
    !print*, 'final slopes 644  :',slope644




    ! why we need these parameter : slope_644, yint644, slope_466, yint466
    ! becasue ==>
    ! Rho_surf_644=slope_644*Rho_surf2p3 + yint644
    ! Then,   ==>
    ! Rho_surf_466=slope_466*Rho_surf644 + yint466

    !********************************************************

    call land_lut % read_land_lut ( path = inp % path)

    call land_lut % sub_land_table (inp % sol,inp%sat,inp%azi,inp%lat, inp%lon, inp%month)

    !	print*,inp%month
    !	print*,'shape of atmospheric path reflectance after sol,sat,azi selected nearest index', shape(land_lut % path_refl_x)

    !**********************

    INT_NL=land_lut % path_refl_x
    Fd_NL=land_lut % Tdn_x
    T_NL=land_lut % Tup_x
    SBAR_NL=land_lut % Sbar_x

    OPTH_NL=land_lut % opt_land_x ! will use / modidy later

    !*********************
    ! land surface elevation correction

    call INT_ELEV(EQWAV_NL,INT_NL,Fd_NL,T_NL,OPTH_NL, &
    SBAR_NL,REF_RAY_NL,MHGHT)



    !***********************************************************************

    !  allocate with n_opt
    allocate ( refl_toa (  land_lut%N_OPT, N_FMR, N_BANDS))

    allocate ( refl_reference ( land_lut% N_OPT))
    n_opt = land_lut % n_opt


    ! Apparent reflectance is calculated from formula
    !
    !  toa_refl = path + (Td Tu Rho/ (1-S Rho))
    !
    !

    allocate ( RHO_surf_2p3(N_OPT,2))


    do i_opt=1, N_OPT


      do i_size = 1,2

        !****************************
        ! estimate surface reflectance here
        ! toa_refl = path + (Fdn Tup Rho/ (1-S Rho))
        ! toa_refl - path  = (Fdn Tup Rho/ (1-S Rho))
        !
        ! P1=toa_refl - path
        ! P2= Tup Fdn
        !
        ! P1 (1-S Rho) = (P2 Rho)
        ! P1 - P1 S Rho = (P2 Rho)
        ! P1 = P1 S Rho + P2 Rho
        ! P1 = Rho ( P1 S + P2)
        ! Rho = P1 /(P1 S + P2)


        P1=inp % rfl(6) -INT_NL(i_opt,i_band2p3,i_size)

        P2=T_NL(i_opt,i_band2p3,i_size)*Fd_NL(i_opt,i_band2p3,i_size)

        RHO_surf_2p3(i_opt,i_size)=P1/((P1*SBAR_NL(i_opt,i_band2p3,i_size))+P2);

        if (RHO_surf_2p3(i_opt,i_size).le.0.01)then
          RHO_surf_2p3(i_opt,i_size)=0.01;
        end if





        ! bright surface
        if(inp % rfl(6).gt.0.18)then
          slope644=slope644+0.025
          slope466=slope466+0.025
        end if



        ! Estiamte surface reflectance using Red/IR and Blue/Red surface correlation
        RHO_surf_644=slope644*RHO_surf_2p3(i_opt,i_size)+yint644;
        RHO_surf_466=slope466*RHO_surf_644+yint466;

        !	print*,RHO_surf_466
        !	print*,RHO_surf_644

        ! Compute model differentiated toa reflectance for band 1 (index 1)
        numb=Fd_NL(i_opt,1,i_size)*T_NL(i_opt,1,i_size)*RHO_surf_466;
        denb=1-SBAR_NL(i_opt,1,i_size)*RHO_surf_466;
        lut_toa_ref(i_opt,1,i_size)=INT_NL(i_opt,1,i_size)+numb/denb; ! LUT

        ! Compute model differentiated toa reflectance for band 3 (index 2)
        numb=Fd_NL(i_opt,2,i_size)*T_NL(i_opt,2,i_size)*RHO_surf_644;
        denb=1-SBAR_NL(i_opt,2,i_size)*RHO_surf_644;
        lut_toa_ref(i_opt,2,i_size)=INT_NL(i_opt,2,i_size)+numb/denb; ! LUT

        ! Compute model differentiated toa reflectance for band 6 (index 3)
        numb=Fd_NL(i_opt,3,i_size)*T_NL(i_opt,3,i_size)*RHO_surf_2p3(i_opt,i_size);
        denb=1-SBAR_NL(i_opt,3,i_size)*RHO_surf_2p3(i_opt,i_size);
        lut_toa_ref(i_opt,3,i_size)=INT_NL(i_opt,3,i_size)+numb/denb; ! LUT


      end do
    end do



    ! loop over fine mode ratio
    do i_fmr = 1,N_FMR
      fmr = (i_fmr-1)/10.

      do i_opt = 1, N_OPT

        do i_cha = 1,n_BANDS
          refl_toa (i_opt,i_fmr,i_cha) = fmr * lut_toa_ref(i_opt,i_cha,1) &
          & + (1 - fmr) * lut_toa_ref(i_opt,i_cha,2)




        end do
      end do




      !- reference channel is #1
      refl_reference = refl_toa(:,i_fmr,CHANNEL_REFERENCE)

      !print*,'refl_reference',refl_reference

      !print*,'input band 1',inp % rfl(CHANNEL_REFERENCE)
      !print*,' AOT array', land_lut % aot_550nm


      !! interpolated in log scale

      aot_temp(i_fmr) = interp1d(refl_reference &
      , land_lut % aot_550nm &
      , inp % rfl(CHANNEL_REFERENCE) &
      , bounds_error = .false. &
      , FILL_VALUE = -999.)


      ! print*,'aot',aot_temp(i_fmr)




    end do


    !print*,'aot_temp',aot_temp



    do i_cha = 1,3
      do i_fmr =1,11

        refl_corrsp(i_fmr,i_cha) = interp1d(land_lut % aot_550nm,refl_toa(:,i_fmr,i_cha) &
        , aot_temp(i_fmr),bounds_error = .false., FILL_VALUE = -999. )

        !  print*,refl_toa(:,i_fmr,i_cha), aot_temp(i_fmr)
        !  print*,'refl_corrsp: ',i_cha,i_fmr, refl_corrsp(i_fmr,i_cha), inp % rfl(i_cha)



      end do

    end do


    do i_cha = 1,3

      n1 = inp % rfl(band(i_cha)) - refl_corrsp(:,i_cha)
      n2 = inp % rfl(band(i_cha))


      err(:,i_cha) = n1/n2
      ! print*,'=====>',i_cha,inp % rfl(band(i_cha)),refl_corrsp(:,i_cha)



    end do


    !  Determine error at 644nm
    !   err_sqrt to err_640 since over land use differntely

    err_640 = err(:,2)



    ! print*,'error :', err_640

    val= minval(abs(err_640))

    idx=minloc(abs(err_640))

    !print*,'idx', idx
    ! print*,'aot',aot_temp(idx(1))
    aod_nleta_temp = aot_temp(idx(1))
    err_nleta_temp = val
    fmf_nleta_temp = (idx(1) - 1 ) /10.
    ! print*,i_fm,i_cm,idx(1)
    aod_allbands(:) = aot_allbands(idx(1),:)




    !val2 = minval(err_nleta_temp)
    !idx2 = minloc(err_nleta_temp)



    !- once find a good match between fwd and measurment give output
    out% aot = aod_nleta_temp
    out % fm_mode =  land_lut % land_fine_mode
    out % cm_mode = 5
    ! print*,fmf_nleta_temp
    !stop
    out % fmf = fmf_nleta_temp
    out% err_n = val


    opt_0470=fmf_nleta_temp * OPTH_NL(:,1,1)+(1-fmf_nleta_temp) * OPTH_NL(:,1,2)
    opt_0640=fmf_nleta_temp * OPTH_NL(:,2,1)+(1-fmf_nleta_temp) * OPTH_NL(:,2,2)
    !	opt_2300=fmf_nleta_temp * OPTH_NL(:,3,1)+(1-fmf_nleta_temp) * OPTH_NL(:,3,2)



    opt_b1 = interp1d(land_lut % aot_550nm,opt_0470 &
    , aod_nleta_temp,bounds_error = .false., FILL_VALUE = -999. )

    opt_b3 = interp1d(land_lut % aot_550nm,opt_0640 &
     , aod_nleta_temp,bounds_error = .false., FILL_VALUE = -999. )

    out% angstrom_exponent= - (log(opt_b1/opt_b3)/log(wavb1/wavb3))



  end subroutine  muri_land_algorithm

  !!**********************************************************************

  subroutine INT_ELEV(EQWAV_NL,INT_NL,Fd_NL,T_NL,OPTH_NL,&
    SBAR_NL,REF_RAY_NL,MHGHT)

    !DESCRIPTION:  Subroutine INTELEV interpolates the lookup
    !               reflectances to the target elevation.
    !               Basically it is fudge:
    !               Interpolate between wavelengths to simulate
    !               elevation by a longer wavelength

    ! Developed by MODIS Aerosol Retrieval Team at NASA GSFC, Greenbelt, MD
    !
    !   INPUT PARAMETERS:
    !         WAV_NL           wavelengths
    !         INT_NL		  radiance
    !	       Fd_NL		  flux down
    !         T_NL         transmission
    !	       SBAR_NL	     sbar
    !         OPTH_NL      optical depth
    !         MGHTH           elevation
    !
    !!OUTPUT PARAMETERS
    !	       INT_NL		  interpolated radiance
    !	       Fd_NL		          interpolated flux down
    !         T_NL                    interpolated transmission
    !	       FdT_NL		  interpolated transmission
    !	       SBAR_NL	  interpolated sbar
    !         OPTH_NL          interpolated optical depth
    !         REF_RAY        Reflectance for rayleigh only
    !         EQWAV_NL        Interpolated wavelengths

    !  real :: surf_alt_out
    !  real :: MHGT  !  MHGHT   Topographic altitude (km)


    !REVISION HISTORY:
    ! 02/08/2006 Shana Mattoo/Rob Levy
    ! Initial revision
    !
    !! TEAM-UNIQUE HEADER:
    !
    ! Developed by MODIS Aerosol Retrieval Team at NASA GSFC, Greenbelt, MD


    IMPLICIT NONE
    SAVE

    !     Inputs
    REAL  :: MHGHT
    integer, parameter :: NLTAU=8
    integer, parameter :: NLWAV=3
    integer, parameter :: NLSIZE=2

    !     Inputs and Outputs
    REAL :: INT_NL(NLTAU,NLWAV,NLSIZE)
    REAL :: Fd_NL(NLTAU,NLWAV,NLSIZE), T_NL(NLTAU,NLWAV,NLSIZE)
    REAL :: SBAR_NL(NLTAU,NLWAV,NLSIZE) !,OPTH_NL(NLTAU,NLWAV,NLSIZE)
    REAL :: REF_RAY_NL(NLWAV),OPTH_NL(NLTAU,NLWAV,NLSIZE)

    REAL :: INT_NL9(NLTAU,NLWAV,NLSIZE)
    REAL :: Fd_NL9(NLTAU,NLWAV,NLSIZE), T_NL9(NLTAU,NLWAV,NLSIZE)
    REAL :: SBAR_NL9(NLTAU,NLWAV,NLSIZE),OPTH_NL9(NLTAU,NLWAV,NLSIZE)

    REAL :: ROD_1013(NLWAV)
    REAL :: p, expfactor
    REAL,PARAMETER :: p0 = 1013.0
    REAL ROD_PRES(NLWAV)
    REAL EQWAV_NL(NLWAV)

    REAL ::lambda0, lambda1, lambda2, diff0, exp0, exp1, exp2
    REAL :: tau0, tau1, tau2

    !!     Dummy
    INTEGER IJ,ISIZE,ITAU,IWAV,IWAV1,IWAV2,JWAV
    REAL  X(2),Y(2),W(2),V(2),Z(2),T(2),U(2)
    REAL  Y1,W1,V1,Z1,T1,U1
    INTEGER LL
    real :: WAV_NL(3)
    integer :: ii

    WAV_NL=(/0.47, 0.64, 2.3/)



    !!     Estimate surface pressure (hypsometric EQ)
    p = p0 * exp(-(MHGHT/7.5))

    !print*,'surface pressure at that elevation',p,p0,MHGHT

    !!     Estimate ROD at nominal wavelenghts at p0 and at pres

    !   print*,'(IWAV) =1 is 0.47um, 2 is 0.64um, 3 is 2.3um'
    DO IWAV = 1, 3
      expfactor = -4.15 + (0.2 * WAV_NL(IWAV))
      ROD_1013(IWAV) = 0.0088 * WAV_NL(IWAV)**expfactor
      ROD_PRES(IWAV) = 0.0088*(p/p0) * WAV_NL(IWAV)**expfactor
      !	print*,'IWAV is ',IWAV
      !	print*,'ROD_1013(IWAV)',ROD_1013(IWAV)
      !    print*,'ROD_PRES(IWAV) corrected Rayleigh Opitcal Depth',ROD_PRES(IWAV)
      !!    Estimate equivalent wavelengths for ROD at pressure
      lambda1 = 0.1
      lambda2 = 4.0
      diff0 = 99.


      ii = 0


      DO WHILE (diff0 .GT. 0.00001)


        lambda0 = (lambda1 + lambda2) / 2.
        exp0 = -4.15 + 0.2*lambda0
        exp1 = -4.15 + 0.2*lambda1
        exp2 = -4.15 + 0.2*lambda2
        tau0 = 0.0088*lambda0**exp0
        tau1 = 0.0088*lambda1**exp1
        tau2 = 0.0088*lambda2**exp2


        IF (tau1 .GT. ROD_PRES(IWAV) .AND.&
        tau2 .LT. ROD_PRES(IWAV)) THEN
        IF (tau0 .GT. ROD_PRES(IWAV)) THEN
          lambda1 = (lambda1 + lambda2)/2.
        ELSE
          lambda2 = (lambda1 + lambda2)/2.
        END IF
      END IF ! tau1 .GT. ROD_PRES(IWAV) .AND.&
      diff0 = ABS(tau0 - ROD_PRES(IWAV))



      ii = ii + 1

      if ( ii .gt. 250) then
        print*,'stop in WHILE loop of land retrieval'
        print*,'diff0',diff0
        print*,'tau0',tau0
        print*,'ROD_PRES(IWAV)',ROD_PRES(IWAV)
        print*,'lambda1',lambda1
        print*,'lambda2',lambda2
        print*,'MHGHT',MHGHT
        stop 'WRONG'
      end if

    END DO ! while loop
    EQWAV_NL(IWAV) = lambda0
  END DO ! IWAV = 1, 3

  !   print*,' Equivalent wavelength of 0.47, 0,64, 2.3 after elevation correction ',EQWAV_NL

  !     INterpolate lookup tables to equiv Waves (let's start only with
  !      1st two wavelengths until we derive 0.86 table)


  DO IJ = 1,2
    X(IJ)=0.0
    Y(IJ)=0.0
    V(IJ)=0.0
    W(IJ)=0.0
    Z(IJ)=0.0
    U(IJ)=0.0
    T(IJ)=0.0
  END DO
  DO 5 ISIZE= 1,NLSIZE
    DO 15 IWAV=1,3
      DO 25  ITAU  = 1,NLTAU
        LL=0
        Y1=0.0
        W1=0.0
        Z1=0.0
        V1=0.0
        IF (IWAV .EQ. 3) THEN
          IWAV1 = IWAV-1
          IWAV2 = IWAV
        ELSE
          IWAV1 = IWAV
          IWAV2 = IWAV+1
        END IF
        DO 60  JWAV = IWAV1, IWAV2
          LL=LL+1

          !    Interpolate on log log scale

          X(LL)=ALOG(WAV_NL(JWAV))
          Y(LL)=ALOG(INT_NL(ITAU,JWAV,ISIZE))
          !                W(LL)=ALOG(FdT_NL(ITAU,JWAV,ISIZE))
          Z(LL)=ALOG(SBAR_NL(ITAU,JWAV,ISIZE))
          V(LL)=OPTH_NL(ITAU,JWAV,ISIZE)
          U(LL)=ALOG(Fd_NL(ITAU,JWAV,ISIZE))
          T(LL)=ALOG(T_NL(ITAU,JWAV,ISIZE))
          IF (OPTH_NL(ITAU,JWAV,ISIZE) .GT. 0.) THEN
            V(LL)=ALOG(OPTH_NL(ITAU,JWAV,ISIZE))
          ENDIF
          60            CONTINUE

          CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,Y,Y1)
          !             CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,W,W1,1)
          CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,Z,Z1)
          CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,V,V1)
          CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,U,U1)
          CALL INTERP_EXTRAP(LL,ALOG(EQWAV_NL(IWAV)),X,T,T1)

          INT_NL9(ITAU,IWAV,ISIZE) = EXP(Y1)
          !             FdT_NL9(ITAU,IWAV,ISIZE) = EXP(W1)
          Fd_NL9(ITAU,IWAV,ISIZE) = EXP(U1)
          T_NL9(ITAU,IWAV,ISIZE) = EXP(T1)
          SBAR_NL9(ITAU,IWAV,ISIZE) = EXP(Z1)
          OPTH_NL9(ITAU,IWAV,ISIZE) = EXP(V1)
          IF (V1 .EQfp. 0.) THEN
            OPTH_NL9(ITAU,IWAV,ISIZE) = V1
          ENDIF


          25          CONTINUE
          15        CONTINUE
          5      CONTINUE


          DO ISIZE= 1,NLSIZE
            DO IWAV=1,3
              DO ITAU  = 1,NLTAU
                INT_NL(ITAU,IWAV,ISIZE) = INT_NL9(ITAU,IWAV,ISIZE)
                Fd_NL(ITAU,IWAV,ISIZE) = Fd_NL9(ITAU,IWAV,ISIZE)
                T_NL(ITAU,IWAV,ISIZE) = T_NL9(ITAU,IWAV,ISIZE)
                SBAR_NL(ITAU,IWAV,ISIZE) = SBAR_NL9(ITAU,IWAV,ISIZE)
                OPTH_NL(ITAU,IWAV,ISIZE) = OPTH_NL9(ITAU,IWAV,ISIZE)
              END DO
              REF_RAY_NL(IWAV) = INT_NL(1,IWAV,ISIZE)
            END DO
          END DO





        end subroutine INT_ELEV



        !*********************************************************************

        subroutine INTERP_EXTRAP(M,X1,X,Y,Y1)

          !----------------------------------------------------------------------
          ! !F77
          !C
          ! !DESCRIPTION: This subroutine is a general purpose routine and
          !              interpolates linearly.  Value of y1 is interpolated
          !              or extrapolated for x1
          !
          ! !INPUT PARAMETERS:I,M,X1,X,Y
          !
          ! !OUTPUT PARAMETERS:Y1,LOPT
          !
          ! !REVISION HISTORY:
          ! $Log: Process_ocean_V2.f,v $
          !
          !
          ! !TEAM-UNIQUE HEADER:
          !
          ! !END
          !----------------------------------------------------------------------

          IMPLICIT NONE

          INTEGER IL,LL,LOPT,M
          REAL PINTEN,PPHI,SINTEN,SPHI
          REAL  X(M),Y(M),Y1,X1,Diff

          Y1=0.0
          LOPT=0
          LL=M-1
          DO 230 IL=1,LL
            !        Extrapolation on lower bound
            IF(X1 .LE.X(1))THEN
              PPHI=X(1)
              SPHI=X(2)
              PINTEN=Y(1)
              SINTEN=Y(2)
              Diff=(SPHI-PPHI)
              if(Diff.LE. 0.0) Diff=1

              Y1=PINTEN+((SINTEN-PINTEN)*((X1-PPHI)/Diff))
              LOPT=1

              !*/  Modified by JC Guu  01/09/97
              !*/  "GO TO 290" is changed to RETURN

              RETURN

              !        Extrapolation on upper bound
            ELSEIF(X1 .GE.X(LL+1)) THEN
              PPHI=X(LL)
              SPHI=X(LL+1)
              PINTEN=Y(LL)
              SINTEN=Y(LL+1)
              Diff=(SPHI-PPHI)
              if(Diff .Le. 0.0) Diff=1
              Y1=PINTEN+((SINTEN-PINTEN)*((X1-PPHI)/Diff))
              LOPT=1

              !*/  Modified by JC Guu  01/09/97
              !*/  "GO TO 290" is changed to RETURN

              RETURN

              !       interpolation
            ELSEIF (X1 .GE.X(IL) .AND.X1.LE. X(IL+1)) THEN
              PPHI=X(IL)
              SPHI=X(IL+1)
              PINTEN=Y(IL)
              SINTEN=Y(IL+1)
              Diff=(SPHI-PPHI)
              if(Diff .Le. 0.0) Diff=1
              Y1=PINTEN+((SINTEN-PINTEN)*((X1-PPHI)/Diff))
              LOPT=1

              !*/  Modified by JC Guu  01/09/97
              !*/  "GO TO 290" is changed to RETURN
              !*/  and two lines are commented out.

              RETURN
              !        ELSE
              !        GO TO 230

            ENDIF
            230    CONTINUE


          end subroutine INTERP_EXTRAP


          !***********************************************************************************

        end module muri_land_retrieval_mod





        !
