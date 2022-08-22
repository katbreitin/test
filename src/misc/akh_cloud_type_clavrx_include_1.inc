!$Id: akh_cloud_type_clavrx_include_1.inc 3190 2019-03-11 22:04:41Z heidinger $
!---------------------------------------------------------------------
! UNIVERSAL_CLOUD_TYPE include file #1 : CLAVR-x specific variables
!----------------------------------------------------------------------
  integer(kind=INT4), PRIVATE, PARAMETER :: Chan_Idx_16um = 6    !channel number for 1.6 micron
  integer(kind=INT4), PRIVATE, PARAMETER :: Chan_Idx_375um = 20  !channel number for 3.75 micron
  integer(kind=INT4), PRIVATE, PARAMETER :: Chan_Idx_67um = 27   !channel number for 6.7 micron
  integer(kind=INT4), PRIVATE, PARAMETER :: Chan_Idx_73um = 28   !channel number for 7.3 micron
  integer(kind=INT4), PRIVATE, PARAMETER :: Chan_Idx_85um = 29   !channel number for 8.5 micron
  integer(kind=INT4), PRIVATE, PARAMETER :: Chan_Idx_104um = 38  !channel number for 10.4 micron
  integer(kind=INT4), PRIVATE, PARAMETER :: Chan_Idx_11um = 31   !channel number for 11 micron
  integer(kind=INT4), PRIVATE, PARAMETER :: Chan_Idx_12um = 32   !channel number for 12 micron
  integer(kind=INT4), PRIVATE, PARAMETER :: Chan_Idx_133um = 33  !channel number for 13.3 micron
  integer(kind=INT4), PRIVATE, SAVE:: Chan_Idx_IRWIN   !channel number for ir window
