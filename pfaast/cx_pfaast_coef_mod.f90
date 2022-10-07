!>
!   $Id: cx_pfaast_coef_mod.f90 3599 2019-11-20 18:49:25Z yli $
!!  @author Andi Walther
!! 
module cx_pfaast_coef_mod
   implicit none
   type coef_type
      real, allocatable :: dry(:,:,:)
      real, allocatable :: ozon(:,:,:)
      real, allocatable :: wvp_cont(:,:,:)
      real, allocatable :: wvp_liquid(:,:,:)
      real, allocatable :: wvp_solid(:,:,:)
      character(len=50) :: sensor  = 'not_set' 
      integer , allocatable :: modis_channel_eqv (:)
      integer , allocatable :: native_channel (:)
      character ( len = 10), allocatable :: channel_descr (:)
   contains
   
      procedure :: allocate_it
      procedure :: deallocate_it
      procedure :: read_it => read_coef 
      procedure :: read_general 
      procedure :: read_modis
      procedure :: read_ahi_abi
      procedure :: read_agri
      procedure :: read_avhrr
      procedure :: read_mersi
      procedure :: read_data_standard
   end type coef_type 
   
   ! - parameter dimensions in coeff files
   integer , parameter :: NK =  5
   integer , parameter :: NL = 101
   integer , parameter :: NM = NL - 1
  
   
   integer , parameter :: NXD = 8
   integer , parameter :: NCD = NXD + 1
   integer , parameter :: LENCD = NCD * NM
   integer , parameter :: LENCDB = LENCD * 4
   
   integer , parameter :: NXO = 9
   integer , parameter :: NCO = NXO + 1
   integer , parameter :: LENCO = NCO * NM
   integer , parameter :: LENCOB = LENCO * 4
   
   integer , parameter :: NXC = 4
   integer , parameter :: NCC = NXC + 1
   integer , parameter :: LENCC = NCC * NM
   integer , parameter :: LENCCB = LENCC * 4   
   
   integer , parameter :: NXL = 2
   integer , parameter :: NCL = NXL + 1
   integer , parameter :: LENCL = NCL * NM
   integer , parameter :: LENCLB = LENCL * 4      
   
   integer , parameter :: NXS = 11
   integer , parameter :: NCS = NXS + 1
   integer , parameter :: LENCS = NCS * NM
   integer , parameter :: LENCSB = LENCS * 4 
   
   integer , parameter :: NXW = NXL + NXS    
   
   character ( len = 3 ) ,parameter ::  comp(NK) = ['dry','ozo','wco','wtl','wts']
   
   character ( len =1020) :: pfast_path
   
   integer, parameter, dimension(5) :: LENCF = [ lencdb,lencob,lenccb,lenclb,lencsb ]
   
   
contains

   !>   Read_coef retruens the coefficients for Pfaast RTM
   !!
   !!
   subroutine read_coef ( this , sensor, ancil_path)
      use cx_string_tools_mod,only: &
         split
      
      implicit none   
      class ( coef_type ) :: this
      character ( len =*) :: sensor
      character ( len =*) :: ancil_path
      character ( len = 15) :: satellite
      character ( len= 15 ) :: device
 
      if ( trim(sensor) == trim(this % sensor)  ) return
      this % sensor = sensor
      
!      print*,'start reading coef ', sensor 
     
      if ( index(sensor,'-') /= 0 ) then
         satellite=sensor
         call split ( satellite,'-', device )  
         
      else
         device = sensor
         satellite = sensor
      end if
      
      pfast_path = trim( ancil_path)//"/static/pfaast/"

      select case ( trim(device) )
      case ('COMS')
         call this % read_general( device, satellite)
      case ('MODIS')
         call this % read_modis( satellite)
      case ('AHI8')
         call this % read_ahi_abi (device, satellite)  
      case ('AHI9')
         call this % read_ahi_abi (device, satellite)
      case ('FY3')
         call this % read_mersi ( device, satellite) 
      case ('FY4')
         call this % read_agri (device, satellite)
      case ( 'AVHRR','HIRS')
         call this % read_avhrr ( device,  satellite )
      case ( 'SEVIRI','MTSAT','VIIRS','FY2','GOES')
         ! - GOES-16-19 is ABI 
         if ( trim(device) .EQ. 'GOES' .and. &
              (satellite .EQ. '16' .OR. satellite .EQ. '17' .OR. satellite .EQ. '18')) then
            call this % read_ahi_abi(device, satellite)
         else 
            call this % read_general ( device, satellite  )
         end if
      case default
         print*,'sensor not found in cx_pfaast:  ',sensor,'with device number/feature: ',device   
         print*,'stopping ...'
         stop    
      end select

   end subroutine read_coef

   !>
   !!
   !!
   subroutine allocate_it (this , n_chan)
      class ( coef_type) :: this
      integer, intent(in) :: n_chan
      
      call this % deallocate_it ()
      
      allocate ( this % dry (NCD, NM, n_chan))
      allocate ( this % ozon (NCO, NM, n_chan))
      allocate ( this % wvp_cont ( NCC, NM , n_chan))
      allocate ( this % wvp_liquid ( NCL, NM , n_chan))
      allocate ( this % wvp_solid ( NCS, NM , n_chan))
      allocate ( this % modis_channel_eqv ( n_chan))
      allocate ( this % native_channel ( n_chan))
      allocate ( this % channel_descr ( n_chan))
      this % native_channel = -1
      this % modis_channel_eqv = -1
   end subroutine allocate_it
   
   !>
   !!
   !!
   subroutine deallocate_it ( this )
      class ( coef_type) :: this
      
      if ( allocated (this % dry) ) deallocate ( this % dry )
      if ( allocated (this % ozon) ) deallocate ( this % ozon )
      if ( allocated (this % wvp_cont) ) deallocate ( this % wvp_cont )
      if ( allocated (this % wvp_liquid) ) deallocate ( this % wvp_liquid )
      if ( allocated (this % wvp_solid) ) deallocate ( this % wvp_solid )
      if ( allocated (this % modis_channel_eqv) ) deallocate ( this % modis_channel_eqv )
      if ( allocated ( this % native_channel) ) deallocate ( this % native_channel )
      if ( allocated (this % channel_descr ) ) deallocate ( this % channel_descr )

      
   
   end subroutine deallocate_it
   
   !> 
   !!  @todo check and test
   !!  @todo make lun safe
   !!
   subroutine read_ahi_abi ( this , sat, sat_num)
      class (coef_type), intent(inout) :: this
      character(len=15), intent(in) :: sat
      character(len=15), intent(in) :: sat_num
      character ( len =100) :: cfile 
      
      integer :: l
      integer, parameter :: ND = 10
      integer, parameter  :: N_CHANNELS = 10   !<  number of channels modis
      real :: coefc ( NCC, NM , N_CHANNELS )
      real :: coefo  (NCO,NM, N_CHANNELS)
      real :: coefd  (ncd,nm, N_CHANNELS)
      real :: coefs  (ncs,nm, N_CHANNELS)
      real :: coefl  (ncl,nm, N_CHANNELS)
    
      integer :: lun_s (NK)
      integer :: krec
      integer :: i,j,k
      
      !- executable

      select case ( trim(sat) )
        case('GOES')
           !--- Will need more logic when goes-18/19 are online.
           if (sat_num .eq. '16') then
             cfile = 'abixxx101.dat'
           else if (sat_num .eq. '17') then
             cfile = 'abi7xxx101.dat'
           else
             cfile = 'abi8xxx101.dat'
           endif
        case('AHI9')
           cfile = 'ahi9xxx101.dat'
        case('AHI8')
          cfile = 'ahixxx101.dat'
        case default
          print*,' wrong name in PAAST .. ',trim(sat) ,' stopping'
          stop
      end select
      
      
      if (sat .eq. 'AHI9' ) then
        call open_files ( cfile, 5, lun_s )
      else if (sat .eq. 'GOES' .and. sat_num .eq. '17') then
        call open_files ( cfile, 5, lun_s )
      else if (sat .eq. 'GOES' .and. sat_num .eq. '18') then
        call open_files ( cfile, 5, lun_s )
      else
        call open_files ( cfile, 4, lun_s )
      endif
      
      ! call open_files ( cfile, 4, lun_s ) 

      do k = 1 , N_CHANNELS
            krec = k 
            read(lun_s(1),rec=krec) ((coefd(i,j,k),i=1,ncd),j=1,nm)
            read(lun_s(2),rec=krec) ((coefo(i,j,k),i=1,nco),j=1,nm)
            read(lun_s(3),rec=krec) ((coefc(i,j,k),i=1,ncc),j=1,nm)
            read(lun_s(4),rec=krec) ((coefl(i,j,k),i=1,ncl),j=1,nm)
            read(lun_s(5),rec=krec) ((coefs(i,j,k),i=1,ncs),j=1,nm)
      end do
         
      do l=1,nk
         close(lun_s(l))
      end do
         
    !  if ( big_endian) then
    !     call flip_rtc(coefd,ncd,nm, N_CHANNELS)
    !     call flip_rtc(coefo,nco,nm, N_CHANNELS)
    !     call flip_rtc(coefc,ncc,nm, N_CHANNELS)
    !     call flip_rtc(coefl,ncl,nm, N_CHANNELS)
    !     call flip_rtc(coefs,ncs,nm, N_CHANNELS)
    !  end if
         
      call this % allocate_it(N_CHANNELS)
         
      ! - return only detector averaged values
      this % dry        = coefd 
      this % ozon       = coefo       
      this % wvp_cont   = coefc 
      this % wvp_solid  = coefs 
      this % wvp_liquid = coefl 
         
      
      
      this % modis_channel_eqv = [ 20,37,27,28,29,30,38,31,32,33 ]
      this % native_channel = [(i , i=7 , 16 ) , 1 ]
        
   end subroutine read_ahi_abi
 
  !  HISTORY:
  !      June 2019: Created by Yue
  !     19 June 2019 (AW): Fixed channel number from 6 to 7  
  !                  and add fake channel 19
  !             At this moment I don't know if channel 19 is correct
  !              Channel 19 is not used by CLAVR-x.
  !
  !
   subroutine read_mersi ( this , sat, sat_num)
      class (coef_type), intent(inout) :: this
      character(len=15), intent(in) :: sat
      character(len=15), intent(in) :: sat_num
      character ( len =100) :: cfile

      integer :: l
      integer, parameter :: ND = 7
      integer, parameter  :: N_CHANNELS = 7   !<  number of channels 
      real :: coefc ( NCC, NM , N_CHANNELS )
      real :: coefo  (NCO,NM, N_CHANNELS)
      real :: coefd  (ncd,nm, N_CHANNELS)
      real :: coefs  (ncs,nm, N_CHANNELS)
      real :: coefl  (ncl,nm, N_CHANNELS)

      integer :: lun_s (NK)
      integer :: krec
      integer :: i,j,k

      !- executable

      if (trim(sat) .ne. 'FY3' .and. sat_num .ne. 'D') then
        print*,' wrong name in PAAST .. ',trim(sat) ,' stopping'
          stop
      end if
      cfile = 'fy3dxxx101.dat'
      call open_files ( cfile, 5, lun_s )

      do k = 1 , N_CHANNELS
            krec = k
            read(lun_s(1),rec=krec) ((coefd(i,j,k),i=1,ncd),j=1,nm)
            read(lun_s(2),rec=krec) ((coefo(i,j,k),i=1,nco),j=1,nm)
            read(lun_s(3),rec=krec) ((coefc(i,j,k),i=1,ncc),j=1,nm)
            read(lun_s(4),rec=krec) ((coefl(i,j,k),i=1,ncl),j=1,nm)
            read(lun_s(5),rec=krec) ((coefs(i,j,k),i=1,ncs),j=1,nm)
      end do

      do l=1,nk
         close(lun_s(l))
      end do


      call this % allocate_it(N_CHANNELS)

      ! - return only detector averaged values
      this % dry        = coefd
      this % ozon       = coefo
      this % wvp_cont   = coefc
      this % wvp_solid  = coefs
      this % wvp_liquid = coefl

      this % modis_channel_eqv = [ 19,20,23,28,29,31,32 ]
      this % native_channel = [(i , i=19 , 25 ) , 1 ]

   end subroutine read_mersi
 
   
   
   !>
   !!  @todo check and test
   !!  @todo make lun safe
   !!
   subroutine read_avhrr ( this , sensor,  satellite )
      class (coef_type), intent(inout) :: this
      character ( len =*) :: sensor
      character ( len =*) :: satellite
      
      character ( len =100) :: cfile = 'avhrncom.dat'      
      
      integer  :: N_CHANNELS   

      integer :: lun_s (NK)
      integer :: i
      integer :: i_sat
      integer :: koff
   
      integer ,parameter :: n_sat = 18
      character ( len =10) , dimension (0 : n_sat -1 ) :: &
             sat_list = ['TIROSN','NOAA06','NOAA07','NOAA08','NOAA09', &
                           'NOAA10','NOAA11','NOAA12','NOAA13','NOAA14', &
                           'NOAA15','NOAA16','NOAA17','NOAA18', &
                           'METOPA','NOAA19','METOPB','METOPC']
   
      ! - Find satellite index from list
      do i_sat = 0 , n_sat -1
        
         if ( trim(satellite) == trim(sat_list ( i_sat )) ) exit
         if ( i_sat == (n_sat -1 )) then
            print*,'PFAAST> Satellite for AVHRR not found! ', trim(satellite)
            print*,'PFAAST> choose one of : ', sat_list
            print*,'PFAAST> stopping'
            stop
         end if   
      end do
      
      ! - WARNING THIS FAKES METOP-C because we don't have it in 
      ! - COEF FILES
      ! - ADDED on April 11 2019 AW
      ! - THIS should be removed ASAP
      
      if (trim(satellite) .eq. 'METOPC') then
        i_sat = i_sat - 1
        print*
        print*,' WARNING: You run METOP-C. We dont have PFAAST coefficients yet.'
        print*,'            This processing uses METOP-B coefficients' 
        print*,'       11 April 2019 (AW)              '
        print*
        
      end if
      !
      
     
      select case ( trim(sensor) )
      
      case ('HIRS')
         n_channels = 19
         koff = i_sat * N_CHANNELS
         cfile = 'hirsccom.dat'
           
      case ('AVHRR')
         n_channels = 3
         cfile = 'avhrncom.dat'
         koff = i_sat * (N_CHANNELS + 1) + 1
         
      
      case default
         print*,'=========== something wron g in pfaast avhrr hirs read ======='
      
      end select
      

      call open_files ( cfile, 6, lun_s ) 
      call this % read_data_standard (n_channels, koff,lun_s)
      
      
            
      select case ( trim(sensor) )
      
      case ('HIRS')
         
!        this % modis_channel_eqv = [ -1,-1,-1,-1,-1,34,33,-1,30,-1,28,27,-1,-1,25,24,-1,23,21 ]
!                                      1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19 ]
!akh mod
         this % modis_channel_eqv = [ -1,-1,-1,36,35,34,33,-1,30,-1,28,27,-1,-1,25,24,-1,23,21 ]
         this % native_channel = [ (i,i=1,19),1 ]   

      case ('AVHRR')
        
         this % modis_channel_eqv = [ 20,31,32 ]
         this % native_channel = [ 3,4,5 ]
      
      case default
         print*,'=========== something wron g in pfaast avhrr hirs read ======='
      
      end select
      
   end subroutine read_avhrr
         
   ! ----------------------------------
   !> This function read MODIS PFAAST coeffecient file
   !!  @param satellite This should be either "AQUA" or "TERRA"
   !!  @todo check channel 26 issue
   !!
   subroutine read_modis ( this , satellite)
         class (coef_type), intent(inout) :: this
     
      character ( len =*) :: satellite
      
      character ( len =100) :: cfile = 'modisdet.com.101.xxx_end'
      
      integer :: l , m
      integer, parameter :: ND = 10
      integer, parameter  :: N_CHANNELS = 17   !<  number of channels modis
      real :: coefc ( NCC, NM , 0:ND, N_CHANNELS )
      real :: coefo  (NCO,NM, 0:ND,N_CHANNELS)
      real :: coefd  (ncd,nm, 0:ND,N_CHANNELS)
      real :: coefs  (ncs,nm, 0:ND,N_CHANNELS)
      real :: coefl  (ncl,nm, 0:ND,N_CHANNELS)
     
      integer :: lun_s (5)
      integer :: krec
      integer :: i,j,k
      integer :: ksat, nsat
      integer, parameter :: NDT = ND + 1 
      integer , parameter :: NRPS = N_CHANNELS * NDT
      integer :: ikrec, krecx
      real :: bufs(lencs)
      
         
      ! define and open the coefficient files
      cfile (18:20) = 'lit'
      !if ( big_endian()) cfile (18:20) = 'big'
      
      ksat = 1 ! terra
      if ( trim(satellite) == 'AQUA' ) ksat = 2
      
      call open_files ( cfile, 10, lun_s )  
 
      
      !  first read cband 26
      !
      !
      ikrec = NRPS * (ksat-1)
      krecx = ikrec + 7
      
      do l = 1, nk
         
         read ( lun_s ( l ) , rec = krecx ) ( bufs(j), j=1, 1)
         nsat = bufs(1)
         if ( nsat /= ksat) then
            print*,'PFAAST> something wrong with MODIS files'
            stop
         end if
         
      
      end do
      
      ! * read in coefficients
      
      krec = ikrec
      do m=0,ND
         do k=1,N_CHANNELS
            krec = krec + 1
            read(lun_s(1),rec=krec) ((coefd(i,j,m,k),i=1,ncd),j=1,nm)
            read(lun_s(2),rec=krec) ((coefo(i,j,m,k),i=1,nco),j=1,nm)
            read(lun_s(3),rec=krec) ((coefc(i,j,m,k),i=1,ncc),j=1,nm)
            read(lun_s(4),rec=krec) ((coefl(i,j,m,k),i=1,ncl),j=1,nm)
            read(lun_s(5),rec=krec) ((coefs(i,j,m,k),i=1,ncs),j=1,nm)
         end do
      end do
      do l=1,nk
         close(lun_s(l))
      enddo
      call this % allocate_it(N_CHANNELS)
      
      ! - return only detector averaged values
      this % dry        = coefd (:,:,0,:)
      this % ozon       = coefo (:,:,0,:)      
      this % wvp_cont   = coefc (:,:,0,:)
      this % wvp_solid  = coefs (:,:,0,:)
      this % wvp_liquid = coefl (:,:,0,:)
      
      this % modis_channel_eqv = [ (i,i=20,36),1 ]
      this % native_channel = [ (i,i=20,36),1 ]
        
   end subroutine read_modis
   
   
   subroutine read_general ( this, sensor,satellite )
      implicit none
      class (coef_type), intent(inout) :: this
      character ( len =*) :: sensor
      character ( len =*) :: satellite
      
      character ( len = 100 ) :: cfile 
      integer :: n_channels 
      integer :: koff
      integer :: lun_s (NK)
      integer :: pos_in_coeffile
      integer :: idx_sat
      integer :: i
      
      select case ( trim ( sensor ))
      
      case ('VIIRS')
          koff =0
          cfile = 'viirscom.dat'
          n_channels = 5
          pos_in_coeffile = 6
          if ( satellite .eq. 'N20') then
             cfile = 'vii20com.dat' 
              n_channels = 7
         end if 
      case ('SEVIRI')
         cfile = 'metsecgencom.dat' 
         n_channels = 8 
         pos_in_coeffile = 10
         idx_sat = 8
         
         if ( satellite == 'MSG09') idx_sat = 9
         if ( satellite == 'MSG10') idx_sat = 10
         if ( satellite == 'MSG11') idx_sat = 11
         koff = (idx_sat - 8 ) * (N_CHANNELS + 1)  + 1

       case ( 'MTSAT')
          cfile = 'mtsatccc_101.dat' 
          N_CHANNELS = 5
          pos_in_coeffile = 6
          read (satellite, * ) idx_sat
          koff = (idx_sat - 1 ) * N_CHANNELS
       
       case('COMS')
          cfile = 'comsccc101.dat'
          N_CHANNELS = 5
          pos_in_coeffile = 5
          koff = 0 
          
       case ('FY2')
         read (satellite, * ) idx_sat
         select case (idx_sat)
         case(1) 
            cfile = 'comsccc101.dat'
         case(2) 
            cfile = 'fy2dccc101.dat'
         case(3) 
            cfile = 'fy2eccc101.dat'
         end select 
         
         N_CHANNELS = 5
         pos_in_coeffile = 5
         koff = 0 
         
      case ('GOES')
         cfile = 'goesccom.dat'
         N_CHANNELS = 23
         pos_in_coeffile = 6
         read (satellite, * ) idx_sat
         koff = (idx_sat - 8 ) * N_CHANNELS
         
      end select
     
      call open_files ( cfile, pos_in_coeffile, lun_s )
     
      call this % read_data_standard (n_channels, koff,lun_s)
     
      select case ( trim ( sensor ))
      case ('VIIRS')
         this % modis_channel_eqv = [20,22,29,31,32]
         this % native_channel = [12,13,14,15,16]
         if ( satellite .eq. 'N20') then
            this % modis_channel_eqv = [-1,-1,20,22,29,31,32]
            this % native_channel = [0,0,12,13,14,15,16]
         end if
      case ('SEVIRI')   
         this % modis_channel_eqv = [20,37,28,29,30,31,32,33]
         this % native_channel = [ (i,i=4,11),1 ]
         
      case ('MTSAT')
         this % modis_channel_eqv = [33,31,32,27,20]     
         this % native_channel = [ 1,2,3,4,5  ]
      case ('FY2')
         this % modis_channel_eqv = [-1,31,32,27,20]      
         this % native_channel = [ 1,2,3,4,5  ] 
      case ('COMS')
         this % modis_channel_eqv = [-1,20,27,31,32]      
         this % native_channel = [ 1,2,3,4,5  ]    
      case ('GOES')
         this % modis_channel_eqv(19:23) = [20,27,31,32,33]
         this % native_channel(19:23) = [ 2 ,3,4,5,6  ] 
      end select
      
      
   end subroutine read_general
   
   

      
   !>   This opens all 5 files
   !!
   !!
   subroutine open_files (cfile, pos , lun_s)
      character ( len =*), intent(in) :: cfile
      integer :: pos
      integer :: lun_s(5)
      character (len =100) :: cfile_loc
      integer :: iux , l
         
      cfile_loc = cfile
      iux = 70

      do l = 1, NK
         cfile_loc(pos : (pos+2)) =comp(l)

         iux = iux + 1
            
         open(iux,file=trim(pfast_path)//trim(cfile_loc),recl=lencf(l), &
                  access='direct', status='old')
         
         lun_s(l)=iux
   
         !--- Output to screen
         !print*,"Opening PFAAST file : ", trim(pfast_path)//trim(cfile_loc)

      end do
      
   end subroutine open_files
      
      !>  THis function reads PFAAST coefficient files for most sensors
      !!
      !!
      subroutine read_data_standard (this, n_channels, koff,  lun_s)
         class (coef_type), intent(inout) :: this
         integer :: n_channels
         integer :: koff
         integer :: lun_s(5)
         integer :: i,j,k,l
         integer :: krec
         
         call this % allocate_it(N_CHANNELS)
         
         
         do k=1,N_CHANNELS
            krec = k + koff
            read(lun_s(1),rec=krec) ((this % dry(i,j,k),i=1,ncd),j=1,nm)   
            read(lun_s(2),rec=krec) ((this % ozon(i,j,k),i=1,nco),j=1,nm)
            read(lun_s(3),rec=krec) ((this % wvp_cont(i,j,k),i=1,ncc),j=1,nm)
            read(lun_s(4),rec=krec) ((this % wvp_liquid(i,j,k),i=1,ncl),j=1,nm)
            read(lun_s(5),rec=krec) ((this % wvp_solid(i,j,k),i=1,ncs),j=1,nm)
         end do
         
          do l=1,nk
            close(lun_s(l))
         enddo
      
      end subroutine read_data_standard

   subroutine read_agri ( this , sat, sat_num)
      class (coef_type), intent(inout) :: this
      character(len=15), intent(in) :: sat
      character(len=15), intent(in) :: sat_num
      character ( len =100) :: cfile 
      
      integer :: l
      integer, parameter :: ND = 8
      integer, parameter  :: N_CHANNELS = 8   !<  number of channels modis
      real :: coefc ( NCC, NM , N_CHANNELS )
      real :: coefo  (NCO,NM, N_CHANNELS)
      real :: coefd  (ncd,nm, N_CHANNELS)
      real :: coefs  (ncs,nm, N_CHANNELS)
      real :: coefl  (ncl,nm, N_CHANNELS)
    
      integer :: lun_s (NK)
      integer :: krec
      integer :: i,j,k
      
      !- executable

      select case ( trim(sat) )
        case('FY4')
           !--- Will need more logic when fy4-b/c/d are online.
           if (sat_num .eq. 'A') then
             cfile = 'fy4xxx101_08.dat'
           endif
        case default
          print*,' wrong name in PFAAST .. ',trim(sat) ,' stopping'
          stop
      end select

      if (sat .eq. 'FY4' .and. sat_num .eq. 'A') then
        call open_files ( cfile, 4, lun_s )
      endif
      
      ! call open_files ( cfile, 4, lun_s ) 

      do k = 1 , N_CHANNELS
            krec = k 
            read(lun_s(1),rec=krec) ((coefd(i,j,k),i=1,ncd),j=1,nm)
            read(lun_s(2),rec=krec) ((coefo(i,j,k),i=1,nco),j=1,nm)
            read(lun_s(3),rec=krec) ((coefc(i,j,k),i=1,ncc),j=1,nm)
            read(lun_s(4),rec=krec) ((coefl(i,j,k),i=1,ncl),j=1,nm)
            read(lun_s(5),rec=krec) ((coefs(i,j,k),i=1,ncs),j=1,nm)
      end do
         
      do l=1,nk
         close(lun_s(l))
      end do
         
    !  if ( big_endian) then
    !     call flip_rtc(coefd,ncd,nm, N_CHANNELS)
    !     call flip_rtc(coefo,nco,nm, N_CHANNELS)
    !     call flip_rtc(coefc,ncc,nm, N_CHANNELS)
    !     call flip_rtc(coefl,ncl,nm, N_CHANNELS)
    !     call flip_rtc(coefs,ncs,nm, N_CHANNELS)
    !  end if
         
      call this % allocate_it(N_CHANNELS)
         
      ! - return only detector averaged values
      this % dry        = coefd 
      this % ozon       = coefo       
      this % wvp_cont   = coefc 
      this % wvp_solid  = coefs 
      this % wvp_liquid = coefl 
         
      this % modis_channel_eqv = [ 33,32,31,29,28,27,20,22 ]
      this % native_channel = [(i , i=1 , 8 ) , 1 ]
        
   end subroutine read_agri
      
      
end module cx_pfaast_coef_mod
