!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: cx_rtm_mod.f90 (src)
!       CX_RTM_MOD (module)
!
!     cx_calculate_rtm ( subroutine)
!
! PURPOSE: This module contains the subroutine which distributes the RTM clear-sky calculations
!
! DESCRIPTION: 
!     aim of this tool is to populate the global variable Trans_Prof_Rtm
!    This is clear-sky transmission profile for IR channels
!
! AUTHORS:
!  Andi Walther, CIMSS, andi.walther@ssec.wisc.edu
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
!--------------------------------------------------------------------------------------

module CX_RTM_MOD

  use NWP_COMMON_MOD,only: &
    nwp_definition
    
  use RTM_COMMON_MOD,only: &
    rtm_params
    
#ifdef LIBRTTOV   
  use CX_RTTOV_BRIDGE_MOD, only: &
    compute_transmission_rttov
#endif

  use CX_PFAAST_MOD, only: &
       COMPUTE_TRANSMISSION_PFAAST

  type cx_rtm_input
    character (len=1024) :: ancil_path
    real ,dimension(:,:),allocatable :: p_std
    real, dimension(:,:), allocatable :: t_prof,w_prof,o_prof, tpw_prof
    real, dimension(:), allocatable ::sat_bin
    integer :: which_rtm
    character (len =20) :: sc_name
    integer :: Chan_Idx
  end type cx_rtm_input
  logical :: first_run = .true.
contains

  subroutine cx_calculate_rtm(inp, trans_prof_rtm)
    implicit none
    type(cx_rtm_input) :: inp
    real,intent(out), dimension(:,:),allocatable :: trans_prof_rtm
    integer :: n_arr(2)
    integer :: ii
    
    n_arr = shape(inp %  p_std)

    if ( .not. (allocated(trans_prof_rtm) ) ) &
        allocate ( trans_prof_rtm(101,n_arr(2)))
  
    trans_prof_rtm = -999.
    if ( inp % which_rtm == 2 ) then
        if (first_run) print*,'Clear Sky Transmission with RTTOV'
#ifdef LIBRTTOV          
      call COMPUTE_TRANSMISSION_RTTOV   ( &
                           trim(inp % Ancil_path) &
                         ,  inp %  p_std &   
                         ,  inp % t_prof &
                         ,  inp % w_prof  &
                         ,  inp % o_prof &
                         ,  inp % Sat_Bin &
                         ,  inp % Sc_Name &
                         ,  inp % Chan_Idx &
                         ,  Trans_Prof_Rtm  &
                         ,  Use_Modis_Channel_Equivalent = .true.  ) 
              
#endif                         
    end if
    
    if (  inp % which_rtm == 1) then
      if (first_run) print*,'Clear Sky Transmission with PFAAST'
      do ii = 1, n_arr(2) 
             
            
         call COMPUTE_TRANSMISSION_PFAAST( &
                           trim(inp % Ancil_path) &
                         ,  inp % t_prof(:,ii) &
                         ,  inp % w_prof(:,ii)  &
                        ,  inp % o_prof(:,ii) &
                         ,  inp % Sat_Bin(ii) &
                        ,  inp % Sc_Name &
                        ,  inp % Chan_Idx &
                        ,  Trans_Prof_Rtm(:,ii) &
                        ,  Use_Modis_Channel_Equivalent = .true.  )
                        
                   
    
       end do
   
    end if
    
    first_run = .false.
    
  
  end subroutine cx_calculate_rtm


end module CX_RTM_MOD
