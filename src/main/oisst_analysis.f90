!$Id: oisst_analysis.f90 3837 2020-05-11 16:05:11Z wstraka $
!--------------------------------------------------------------------------------------
! Clouds from AVHRR Extended (CLAVR-x) 1b PROCESSING SOFTWARE Version 6.0
!
! NAME: oisst_analysis.f90 (src)
!       OISST_ANALYSIS (program)
!
! PURPOSE: Routine to handle the Reynolds OISST analysis
!
! DESCRIPTION: see below
!
! AUTHORS:
!  Andrew Heidinger, Andrew.Heidinger@noaa.gov
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
! REVISION HISTORY:
!  July 2004 - modified to look for previous file if nearest is not present
!  Apr 2007 - Lahey fortran does not all byte swapping at compilation so it requires
!   the CONVERT="BIG_ENDIAN" in the open statement.  This is not standard.
!  October 2009 - moved to daily 0.25 degree Reynolds SST
!
!
! data is one degree resolution
! 
! http://www.emc.ncep.noaa.gov/research/cmb/sst_analysis/
!
! global (89.875S - 89.875N)  (1440 x 720) starts at 89.5S and GM
!
!
!  DESCRIPTION OF THE DAILY OI SEA SURFACE TEMPERATURE (SST) ANALYSIS Version2
!  
!  The SST analysis is computed daily on a 0.25 degree latitude/longitude grid. 
!  This is version 1.0. There are two products with different satellite data. 
!  Both products use in situ data from ships and buoys. Also SSTs are generated 
!  for sea-ice concentrations above 50%. The sea ice for 1981-2004 is from 
!  http://nsidc.org/data/nsidc-0051.html 
!  (Cavalieri D., C. Parkinson, P. Gloerson, and H.J. Zwally. 1997, updated 2005. 
!  Sea ice concentrations from Nimbus-7 SMMR and DMSP SSM/I passive microwave data, 
!  June to September 2001. Boulder, CO, USA). The sea ice from 2005 to present is from 
!  http://polar.ncep.noaa.gov/seaice/ 
!  (Grumbine, R. W., 1996: Automated passive microwave sea ice concentration 
!  analysis at NCEP, 13pp. Unpublished manuscript available from NCEP/NWS/NOAA, 
!  5200 Auth Road, Camp Springs, MD, 20746, USA.)
!  
!  The first product uses NODC's AVHRR Pathfinder Version 5 
!  http://pathfinder.nodc.noaa.gov 
!  for September 1, 1981 though December 31, 2005 and the operational US Navy AVHRR data 
!  (May, D.A., M. M. Parmeter, D. S. Olszewski and B. D. McKenzie, 1998: Operational 
!  processing of satellite sea surface temperature retrievals at the Naval Oceanographic 
!  Office, Bull. Amer. Met. Soc., 79, 397-407) from January 1, 2006, through present. 
!  This product will henceforth be referred to as the AVHRR product.
!  
!  The second product adds AMSR-E version 5 data obtained from 
!  http://www.remss.com/ 
!  along with the AVHRR data used in version 1a and is available from June 1, 2002, 
!  (the start of AMSR-E) through present. This product will henceforth be 
!  termed the AVHRR + AMSR product.
!  
!  Both analyses include a bias correction of the satellite data with respect to 
!  in situ data using an empirical orthogonal teleconnection (EOT) algorithm. 
!  A short description of the complete analysis procedure can be found in the 
!  AMS extended abstract file (Reynolds-reviewed-rev.pdf).
!  
!  The SST analyses are available in individual daily files. The AVHRR product 
!  is named avhrr-only-v2.YYYYMMDD where YYYY is the year, MM is the month, 
!  and DD is the day. The files can be found on 
!  ftp://eclipse.ncdc.noaa.gov/pub/OI-daily-v2/IEEE/YYYY/AVHRR 
!  where YYYY is the year: 1981 to present. The files were written in IEEE 
!  binary (big-endian) and must be decompressed using gunzip. 
!  The AVHRR + AMSR-E product is written with the same format as the AVHRR product.
!  However, the file names are avhrr-only-v2.YYYYMMDD. The files can be found on
!  ftp://eclipse.ncdc.noaa.gov/pub/OI-daily-v2/IEEE/YYYY/AVHRR-AMSR 
!  where YYYY is the year: 2002 to present.
!  
!  Each file contains 4 records with integer*4 year, month, day, 
!  followed by a gridded integer*2 array. The first array is SST. 
!  The second array is the SST anomaly with respect to a 1971-2000 
!  base period. The third array is the sea ice concentration. 
!  The fourth array is the standard deviation of the analysis error 
!  which includes sampling, random and bias error. 
!  
!  Note: The SST, SST ANOMALY AND ERROR ARRAYS MUST BE MULTIPLIED BY 0.01 
!  TO CONVERT THE VALUES TO DEGREE C. The sea ice concentration array is in per cent (0-100). 
!  Missing values are -999.
!  
!  All arrays consist of 1440 spatial points in longitude from 0.125E to 359.875E 
!  in intervals of 0.25 increasing eastward, and 720 spatial points in latitude 
!  from 89.875S to 89.875N in intervals of 0.25 increasing northward.
!  
!  Each day consists of four FORTRAN records:
!  1. Three 4-byte integers for the year, month and day followed by 
!     1440*720 2-byte integer SST values.
!  2. Three 4-byte integers for the year, month and day followed by 
!     1440*720 2-byte integer SST anomaly values.
!  3. Three 4-byte integers for the year, month and day followed by 
!     1440*720 2-byte integer error values.
!  4. Three 4-byte integers for the year, month and day followed by 
!     1440*720 2-byte integer ice concentration values.
!  
!  Each record is written with a FORTRAN unformatted write which adds 
!  an extra 4 byte header and trailer word to the total record.
!--------------------------------------------------------------------------------------
module OISST_ANALYSIS
  use univ_kind_defs_mod, only: i2, f4
  use CONSTANTS_MOD
  use CX_DATE_TIME_TOOLS_MOD
  use PIXEL_COMMON_MOD
  use FILE_UTILS, only: file_test
  use CX_NETCDF4_MOD
  use CLAVRX_MESSAGE_MOD

  implicit none

  private

  public:: GET_PIXEL_SST_ANALYSIS, GET_OISST_MAP_FILENAME, READ_OISST_ANALYSIS_MAP

  integer, parameter, private:: num_lon_sst_anal =1440, num_lat_sst_anal= 720 
  real, parameter, private:: first_lon_sst_anal = 0.125, first_lat_sst_anal = -89.875, & 
                            last_lon_sst_anal = 359.875, last_lat_sst_anal = 89.875
  real, parameter, private:: del_lon_sst_anal = (last_lon_sst_anal - first_lon_sst_anal)/(num_lon_sst_anal-1)
  real, parameter, private:: del_lat_sst_anal = (last_lat_sst_anal - first_lat_sst_anal)/(num_lat_sst_anal-1)
  real, parameter, private:: min_sst_anal = -4.0, max_sst_anal=36.0

  integer(kind=i2), dimension(num_lon_sst_anal,num_lat_sst_anal,1,1) ::  &
       temp_i2_buffer
  real(kind=real4), dimension(num_lon_sst_anal,num_lat_sst_anal), private, save:: oisst_anal_map
  real(kind=real4), dimension(num_lon_sst_anal,num_lat_sst_anal), private, save:: oisst_err_map
  real(kind=real4), dimension(num_lon_sst_anal,num_lat_sst_anal), private, save:: oisst_anal_map_uni
  real(kind=real4), dimension(num_lon_sst_anal,num_lat_sst_anal), private, save:: oisst_cice_map

  integer(kind=int4), parameter, private :: MAX_OISST_LATENCY = 5 !including current day

  contains

  !-------------------------------------------------------------------
  ! Function to find the oisst map name.
  !-------------------------------------------------------------------
  function GET_OISST_MAP_FILENAME(year_in,day_of_year,oisst_path) result(oisst_filename)

   character(*), intent(in) :: oisst_path
   integer(kind=int2), intent(in):: year_in
   integer(kind=int2), intent(in):: day_of_year
  
   character(len=1020) :: oisst_filename
   character(len=1020) :: oisst_filename_tmp
   character(len=1020) :: oisst_filename_tmp_preliminary
   character(len=1020) :: oisst_filename_v21
   character(len=1020) :: oisst_filename_tmp_v21
   character(len=1020) :: oisst_filename_tmp_preliminary_v21
   integer(kind=int4) :: iday, year, month, day, jday, ileap
   character (len=2)   :: day_string, month_string
   character (len=4)   :: year_string

    oisst_filename = "no_file"

   do iday=0, MAX_OISST_LATENCY - 1
      jday = day_of_year - iday
      year = year_in 
      ileap = leap_year_fct(year)
      
      if (jday < 1) then
         year = year - 1
         ileap = leap_year_fct(year)
         jday = (365 + ileap) + jday
      end if 
       
      month = compute_month(jday, ileap)
      day = compute_day(jday, ileap)
      write (year_string,fmt="(I4)") year
      write (month_string, '(I2.2)') month
      write (day_string,   '(I2.2)') day

      
      ! V2 OISST
      oisst_filename_tmp= "avhrr-only-v2."//year_string//month_string//day_string
      oisst_filename_tmp = trim(oisst_path)//trim(year_string)//"/"//trim(oisst_filename_tmp)
      oisst_filename_tmp_preliminary = trim(oisst_filename_tmp)//"_preliminary.nc"
      oisst_filename_tmp = trim(oisst_filename_tmp)//".nc"

      ! V2r1 OISST
      ! this has a slightly different naming convention than the V2 version, so we need to 
      ! build the preliminary and actual names - WCS3
      oisst_filename_tmp_v21= "oisst-avhrr-v02r01."//year_string//month_string//day_string
      oisst_filename_tmp_v21 = trim(oisst_path)//trim(year_string)//"/"//trim(oisst_filename_tmp_v21)
      oisst_filename_tmp_preliminary_v21 = trim(oisst_filename_tmp_v21)//"_preliminary.nc"
      oisst_filename_tmp_v21 = trim(oisst_filename_tmp_v21)//".nc"
      
      ! This tests to see if either file exists
      if ((file_test(trim(oisst_filename_tmp)) ) .OR. &
           (file_test(trim(oisst_filename_tmp_v21)) ) )then
           
           !default is to use the old naming convention, secondary is the "new" name
           
           if (file_test(trim(oisst_filename_tmp)) ) then
                 oisst_filename = oisst_filename_tmp
                 call MESG("Found "//trim(oisst_filename_tmp), &
                            level=Verb_Lev%VERBOSE)
                exit
           endif
           
           if (file_test(trim(oisst_filename_tmp_v21)) ) then
                 oisst_filename = oisst_filename_tmp_v21
                 call MESG("Found "//trim(oisst_filename_tmp_v21),&
                           level=Verb_Lev%VERBOSE)
                exit
           endif
      end if

      
      !--- check for preliminary file (true of recent files)
      !--- same check as above done for OISST v2.1
      
      if ((file_test(trim(oisst_filename_tmp_preliminary)) ) .OR. &
           (file_test(trim(oisst_filename_tmp_preliminary_v21)) ) ) then
                      
           if (file_test(trim(oisst_filename_tmp_preliminary)) ) then
                 oisst_filename = oisst_filename_tmp_preliminary
                 call MESG("Found "//trim(oisst_filename_tmp_preliminary), &
                            level=Verb_Lev%VERBOSE)
                exit
           endif
           
           if (file_test(trim(oisst_filename_tmp_preliminary_v21)) ) then
                 oisst_filename = oisst_filename_tmp_preliminary_v21
                 call MESG("Found "//trim(oisst_filename_tmp_preliminary_v21), &
                            level=Verb_Lev%VERBOSE)
                exit
           endif
      end if

   end do
   
   if (oisst_filename == "no_file") then
      call MESG("WARNING, OISST file read failed", level=Verb_Lev%WARNING)
   endif
      
   return

end function GET_OISST_MAP_FILENAME

!----------------------------------------------------------------------
! Read the Data
!----------------------------------------------------------------------
  subroutine READ_OISST_ANALYSIS_MAP(sst_anal_name)

   character(len=*), intent(in):: sst_anal_name
   integer(kind=int4):: i,j
   integer:: ii,jj
   integer:: ncid
   integer, dimension(4) :: sds_start_4d, sds_edge_4d, sds_stride_4d

   !------------------------------------------------------------------------------
   ! read in data
   !------------------------------------------------------------------------------

   call OPEN_NETCDF(trim(sst_anal_name), ncid)
   if (ncid < 0) then 
     call MESG("Error opening OISST Analysis file, file = "//trim(sst_anal_name), level=Verb_Lev%WARNING)
     use_sst_anal = 0
     return
   endif

   sds_start_4d = 1
   sds_edge_4d(1) = num_lon_sst_anal
   sds_edge_4d(2) = num_lat_sst_anal 
   sds_edge_4d(3) = 1
   sds_edge_4d(4) = 1
   sds_stride_4d = 1

   call READ_NETCDF(ncid, sds_start_4d, sds_stride_4d, sds_edge_4d, 'sst',  &
        temp_i2_buffer)
   oisst_anal_map = real(temp_i2_buffer(:,:,1,1), kind=f4)*0.01_f4

   call READ_NETCDF(ncid, sds_start_4d, sds_stride_4d, sds_edge_4d, 'err',  &
        temp_i2_buffer)
   oisst_err_map = real(temp_i2_buffer(:,:,1,1), kind=f4)*0.01_f4

   call READ_NETCDF(ncid, sds_start_4d, sds_stride_4d, sds_edge_4d, 'ice',  &
        temp_i2_buffer)
   oisst_cice_map = real(temp_i2_buffer(:,:,1,1), kind=f4)*0.01_f4

   call CLOSE_NETCDF(ncid)

   call MESG("OISST analysis read in successfully", level=Verb_Lev%DEFAULT)

   !-------- convert from Celsius to Kelvin
   oisst_anal_map = oisst_anal_map + 273.15

   where(oisst_anal_map < 270.0)
     oisst_anal_map = missing_value_real4
   endwhere

   !---- compute the 3x3 uniformity of sst_anal
   do i = 1, num_lon_sst_anal
    ii = max(2,min(num_lon_sst_anal-1,i))
    do j = 1, num_lat_sst_anal
      jj = max(2,min(num_lat_sst_anal-1,j))
       oisst_anal_map_uni(i,j) = (maxval(oisst_anal_map(ii-1:ii+1,jj-1:jj+1)) - &
                        minval(oisst_anal_map(ii-1:ii+1,jj-1:jj+1)) ) / 3.0
    enddo
   enddo

   where (oisst_anal_map_uni > 8.0)
    oisst_anal_map_uni = 8.0
   endwhere

  end subroutine READ_OISST_ANALYSIS_MAP

!----------------------------------------------------------------------------------------
! Interpolate SST Analysis to each AVHRR pixel
!
!  input: j1 - the first scan index of this segment
!         j2 - the last scan index of this segment
!----------------------------------------------------------------------------------------
subroutine GET_PIXEL_SST_ANALYSIS(j1,j2)
  integer, intent(in):: j1,j2
  integer:: i
  integer:: j
  integer:: ilon_sst_anal
  integer:: ilat_sst_anal
  integer:: ilon_sst_anal_x
  integer:: ilat_sst_anal_x
  real:: xlon
  real:: lon_sst_anal
  real:: lat_sst_anal
  real:: lon_weight
  real:: lat_weight

  !--- initialize output
  sst_anal =  missing_value_real4
  sst_anal_err =  missing_value_real4
  sst_anal_cice =  missing_value_real4
  sst_anal_uni =  missing_value_real4

  do j = j1, j1+j2-1

     do i = 1, Image%Number_Of_Elements

     if (bad_pixel_mask(i,j) == sym%YES) then
      cycle
     endif

     if (nav % lon(i,j) < 0.0) then
         xlon = nav % lon(i,j) + 360.0
     else
         xlon = nav % lon(i,j)
     endif

     ilon_sst_anal = max(1,min(num_lon_sst_anal,int((xlon - first_lon_sst_anal)/del_lon_sst_anal+1)))
     ilat_sst_anal = max(1,min(num_lat_sst_anal,int((nav % lat(i,j) - first_lat_sst_anal)/del_lat_sst_anal+1)))
     ilon_sst_anal_x = ilon_sst_anal+1
     ilat_sst_anal_x = ilat_sst_anal+1

     ilon_sst_anal_x = max(1,min(num_lon_sst_anal,ilon_sst_anal_x))
     ilat_sst_anal_x = max(1,min(num_lat_sst_anal,ilat_sst_anal_x))

     lon_sst_anal = first_lon_sst_anal + (ilon_sst_anal-1) * del_lon_sst_anal
     lat_sst_anal = first_lat_sst_anal + (ilat_sst_anal-1) * del_lat_sst_anal

     lon_weight = max(0.0,min(1.0,(xlon - lon_sst_anal) / del_lon_sst_anal))
     lat_weight = max(0.0,min(1.0,(nav % lat(i,j) - lat_sst_anal) / del_lat_sst_anal))

     if (minval(oisst_anal_map(ilon_sst_anal:ilon_sst_anal_x,ilat_sst_anal:ilat_sst_anal_x))  &
                /= missing_value_real4) then
        sst_anal(i,j) = (1.0-lon_weight)*(1.0-lat_weight)*oisst_anal_map(ilon_sst_anal,ilat_sst_anal) + &
                     (lon_weight)*(1.0-lat_weight)*oisst_anal_map(ilon_sst_anal_x,ilat_sst_anal) + &
                     (1.0-lon_weight)*(lat_weight)*oisst_anal_map(ilon_sst_anal,ilat_sst_anal_x) + &
                     (lon_weight)*(lat_weight)*oisst_anal_map(ilon_sst_anal_x,ilat_sst_anal_x) 
     else

        sst_anal(i,j) = oisst_anal_map(ilon_sst_anal,ilat_sst_anal)

     endif

     sst_anal_err(i,j) =  oisst_err_map(ilon_sst_anal,ilat_sst_anal)
     sst_anal_cice(i,j) = oisst_cice_map(ilon_sst_anal,ilat_sst_anal)
     sst_anal_uni(i,j) = oisst_anal_map_uni(ilon_sst_anal,ilat_sst_anal)

    end do
  end do

  where(sst_anal < 0.0)
    sst_anal = missing_value_real4
  endwhere
  where(sst_anal_cice < 0.0)
    sst_anal_cice = missing_value_real4
  endwhere
  where(sst_anal_uni < 0.0)
    sst_anal_uni = missing_value_real4
  endwhere
  where(sst_anal_err < 0.0)
    sst_anal_err = missing_value_real4
  endwhere

end subroutine GET_PIXEL_SST_ANALYSIS


end module OISST_ANALYSIS
