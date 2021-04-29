!-------------------------------------------------------------------------------
! Get the probabilitz for nb Classifiers
!-------------------------------------------------------------------------------

module NBM_CLOUD_MASK_GET_PROB_MASK_PHASE

 use NB_CLOUD_MASK_SERVICES

 implicit none

 public GET_PROB_MASK_PHASE

 contains

subroutine GET_PROB_MASK_PHASE (X,Y,Z,Satzen, Solzen, Lunzen, Solglintzen, Lunglintzen,&
                        Lunglint_Mask, Moon_Illum_Frac, City_Mask, &
                        Solscatang, Tpw, Tsfc, Zsfc, Zsfc_Std, Land_Class , Snow_Class, Solglint_Mask, &
                        Coast_Mask, Isfc, Class_Idx, Lut, Missing_Value_Int, Missing_Value_Real, &
                        Clear_Cond_Ratio, Water_Cond_Ratio, Ice_Cond_Ratio, Obs_Prob)

! Assumption is that x,y,z are all initialized by call routine.

    real, intent(in) :: X, Y, Z
    real, intent(in) :: Satzen, Solzen, Lunzen
    real, intent(in) :: Solglintzen, Lunglintzen, Moon_Illum_Frac
    real, intent(in) :: Solscatang, Missing_Value_Real, Zsfc, Zsfc_Std
    real, intent(in) :: Tpw, Tsfc
    integer(kind=int1), intent(in) :: Solglint_Mask, Lunglint_Mask, Land_class, Snow_Class, Coast_Mask, City_Mask
    integer(kind=int1),intent(in) :: Missing_Value_Int, Isfc
    integer, intent(in) :: Class_Idx
    real, intent(out) :: Clear_Cond_Ratio, Water_Cond_Ratio, &
                                     Ice_Cond_Ratio, Obs_Prob
    type(Classifier), intent(in), dimension (:) :: Lut
    real :: Obs_Prob_Thresh
    integer :: Ix, Iy, Iz
    integer :: Classifier_Rank
    
    Classifier_Rank = Lut(Class_Idx)%Rank
    
    
    !initialize variables
    Clear_Cond_Ratio = Missing_Value_Real
    Water_Cond_Ratio = Missing_Value_Real
    Ice_Cond_Ratio = Missing_Value_Real
    Obs_Prob = Missing_Value_Real
    
    !-- check for data that is withing the limits for this classifier
    if (Lut(Class_Idx)%On_Flag(Isfc) == 0) return

    if (Lut(Class_Idx)%Zen_Min /= Missing_Value_Real .and. Lut(Class_Idx)%Zen_Max /= Missing_Value_Real) then 
       if ( Satzen /= Missing_Value_Real .and. &
           (Satzen < Lut(Class_Idx)%Zen_Min .OR.Satzen > Lut(Class_Idx)%Zen_Max)) return
    endif

    if (Lut(Class_Idx)%Solzen_Min /= Missing_Value_Real .and. Lut(Class_Idx)%Solzen_Max /= Missing_Value_Real) then 
       if (Solzen /= Missing_Value_Real .and. &
           (Solzen < Lut(Class_Idx)%Solzen_Min .OR. Solzen > Lut(Class_Idx)%Solzen_Max)) return
    endif

    if (Lut(Class_Idx)%Solglintzen_Min /= Missing_Value_Real .and. Lut(Class_Idx)%Solglintzen_Max /= Missing_Value_Real) then 
       if (Solglintzen /= Missing_Value_Real .and. &
           Land_Class >= 5 .and. &
           Snow_Class == 1 .and. &  
           (Solglintzen < Lut(Class_Idx)%Solglintzen_Min .OR. Solglintzen > Lut(Class_Idx)%Solglintzen_Max)) return
    endif

    if (Lut(Class_Idx)%Lunglintzen_Min /= Missing_Value_Real .and. Lut(Class_Idx)%Lunglintzen_Max /= Missing_Value_Real) then 
       if (Lunglintzen /= Missing_Value_Real .and. &
           Land_Class >= 5 .and. &
           Snow_Class == 1 .and. &  
           (Lunglintzen < Lut(Class_Idx)%Lunglintzen_Min .OR. &
            Lunglintzen > Lut(Class_Idx)%Lunglintzen_Max)) return
    endif

    if (Solglint_Mask >= 0 .and. Land_Class >= 5 .and. Snow_Class == 1 .and. &
       (Solglint_Mask < Lut(Class_Idx)%Solglint_Mask_Min .OR. &
        Solglint_Mask > Lut(Class_Idx)%Solglint_Mask_Max)) return

    if (Lunglint_Mask >= 0 .and. Land_Class >= 5 .and. Snow_Class == 1 .and. &
       (Lunglint_Mask < Lut(Class_Idx)%Lunglint_Mask_Min .OR. &
        Lunglint_Mask > Lut(Class_Idx)%Lunglint_Mask_Max)) return


    if (Moon_Illum_Frac /= Missing_Value_Real .and. &
       (Moon_Illum_Frac < Lut(Class_Idx)%Moon_Illum_Frac_Min .or. &
        Moon_Illum_Frac > Lut(Class_Idx)%Moon_Illum_Frac_Max)) return

    if (City_Mask /= Missing_Value_Int .and. &
       (City_Mask < Lut(Class_Idx)%City_Mask_Min .or. &
        City_Mask > Lut(Class_Idx)%City_Mask_Max)) return

    if (Lut(Class_Idx)%Solscatang_Min /= Missing_Value_Real .and. Lut(Class_Idx)%Solscatang_Max /= Missing_Value_Real) then 
       if (Solscatang /= Missing_Value_Real .and. &
          (Solscatang < Lut(Class_Idx)%Solscatang_Min .OR. Solscatang > Lut(Class_Idx)%Solscatang_Max)) return
    endif

    if (Lut(Class_Idx)%Tpw_Min /= Missing_Value_Real .and. Lut(Class_Idx)%Tpw_Max /= Missing_Value_Real) then 
       if (Tpw /= Missing_Value_Real .and. &
          (Tpw < Lut(Class_Idx)%Tpw_Min .OR. Tpw > Lut(Class_Idx)%Tpw_Max)) return
    endif

    if (Lut(Class_Idx)%Tsfc_Min /= Missing_Value_Real .and. Lut(Class_Idx)%Tsfc_Max /= Missing_Value_Real) then 
       if (Tsfc /= Missing_Value_Real .and. &
          (Tsfc < Lut(Class_Idx)%Tsfc_Min .OR. Tsfc > Lut(Class_Idx)%Tsfc_Max)) return
    endif

    if (Lut(Class_Idx)%Zsfc_Min /= Missing_Value_Real .and. Lut(Class_Idx)%Zsfc_Max /= Missing_Value_Real) then 
       if (Zsfc < Lut(Class_Idx)%Zsfc_Min .OR. Zsfc > Lut(Class_Idx)%Zsfc_Max) return
    endif

    if (Lut(Class_Idx)%Zsfc_Std_Min /= Missing_Value_Real .and. Lut(Class_Idx)%Zsfc_Std_Max /= Missing_Value_Real) then 
       if (Zsfc_Std < Lut(Class_Idx)%Zsfc_Std_Min .OR. Zsfc_Std > Lut(Class_Idx)%Zsfc_Std_Max) return
    endif

    if (Snow_Class >= 1 .and. &
       (Snow_Class < Lut(Class_Idx)%Snow_Class_Min .OR. &
        Snow_Class > Lut(Class_Idx)%Snow_Class_Max)) return

    if (Coast_Mask >= 0 .and. &
       (Coast_Mask < Lut(Class_Idx)%Coast_Mask_Min .OR. &
        Coast_Mask > Lut(Class_Idx)%Coast_Mask_Max)) return

    !--- check for valid data
    if (X == Missing_Value_Real) return
    if (Isfc .le. 0) return

    if (Classifier_Rank .ge. 2) then 
       if (Y == Missing_Value_Real) return
    endif

    if (Classifier_Rank .ge. 3) then 
       if (Z == Missing_Value_Real) return
    endif


     
    ! --- determine probabilitz

    Ix = min((Lut(Class_Idx)%Nbins_X), max(1,NINT((X - &
                                           Lut(Class_Idx)%X_Min) / &
                                           Lut(Class_Idx)%X_Bin)))
    Iy = 0
    Iz = 0

    if (Classifier_Rank >=  2) then 
       if (Y == Missing_Value_Real) return
       Iy = min((Lut(Class_Idx)%Nbins_Y), max(1,NINT((Y - &
                                           Lut(Class_Idx)%Y_Min) / &
                                           Lut(Class_Idx)%Y_Bin)))
    endif

    if (Classifier_Rank >= 3) then 
       if (Z == Missing_Value_Real) return
       Iz = min((Lut(Class_Idx)%Nbins_Z), max(1,NINT((Z - &
                                           Lut(Class_Idx)%Z_Min) / &
                                           Lut(Class_Idx)%Z_Bin)))
    endif
    
    Obs_Prob_thresh = 0.0
    
    select case (Classifier_Rank)
        case (1)
            Obs_Prob = Lut(Class_Idx)%Obs_Table(Ix,Isfc,1,1)
            if (Obs_Prob > Obs_Prob_thresh) then
               Clear_Cond_Ratio = Lut(Class_Idx)%Clear_Table(Ix,Isfc,1,1)
               Ice_Cond_Ratio =   Lut(Class_Idx)%Ice_Table(Ix,Isfc,1,1)
               Water_Cond_Ratio = Lut(Class_Idx)%Water_Table(Ix,Isfc,1,1)
            endif
        case (2)
            Obs_Prob = Lut(Class_Idx)%Obs_Table(Ix,Iy,Isfc,1)
            if (Obs_Prob > Obs_Prob_thresh) then
               Clear_Cond_Ratio = Lut(Class_Idx)%Clear_Table(Ix,Iy,Isfc,1)
               Ice_Cond_Ratio =   Lut(Class_Idx)%Ice_Table(Ix,Iy,Isfc,1)
               Water_Cond_Ratio = Lut(Class_Idx)%Water_Table(Ix,Iy,Isfc,1)
            endif
        case (3)
            Obs_Prob = Lut(Class_Idx)%Obs_Table(Ix,Iy,Iz,Isfc)
            if (Obs_Prob > Obs_Prob_thresh) then
               Clear_Cond_Ratio = Lut(Class_Idx)%Clear_Table(Ix,Iy,Iz,Isfc)
               Ice_Cond_Ratio =   Lut(Class_Idx)%Ice_Table(Ix,Iy,Iz,Isfc)
               Water_Cond_Ratio = Lut(Class_Idx)%Water_Table(Ix,Iy,Iz,Isfc)
            endif
        case DEFAULT
            PRINT *, "INVALID CLASSIFIER RANK"
            return
        END select


return

end subroutine GET_PROB_MASK_PHASE

end module NBM_CLOUD_MASK_GET_PROB_MASK_PHASE
