#+ Produces a preliminary APPLICATION-INDEPENDENT configure script that
#  "finishes" a Makefile-format file and a sh-format file
#  ('include_appind.mk.include' and 'include_appind.sh.include') based on
#  templates ('include_appind.mk.include-template' and
#  'include_appind.sh.include-template'). These output files include many
#  useful values, which have been set with the greatest practical script
#  portability in mind.  A C-syntax header file is also produced
#  (cfghdr_appind.h-template -> cfghdr_appind.h), and it may contain cpp
#  "defines", etc.
#
#  Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT([Application-independent information], [1.0], tmichaels@wisc.edu)
AC_REVISION([Revision: 1.0])

# Try to ascertain the general type of this system:
AC_CANONICAL_BUILD

# Check for useful utility programs:
AC_PROG_GREP
AC_PROG_AWK
AC_PROG_SED
AC_PROG_LN_S
AC_PROG_MKDIR_P  # Needs GNU 'install-sh', 'config.sub', and 'config.guess'

# Define some output files, etc.:

AC_CONFIG_FILES([config_appind.mk.include:config_appind.mk.include-template])
AC_CONFIG_FILES([config_appind.sh.include:config_appind.sh.include-template])
AC_CONFIG_HEADERS([cfghdr_appind.h:cfghdr_appind.h-template])

# "Include" build environment configuration file:
CFG_BUILDENV=build_env_config;
if test -f $CFG_BUILDENV; then
   . "./$CFG_BUILDENV";
else
   BLAH=;
fi

failed_aux_libs=""; usable_aux_libs=""; F2003_ability=no;

for proglang in $prog_languages_needed
do
   LDFLAGS=; export LDFLAGS;  LIBS=; export LIBS;
   FCFLAGS=; export FCFLAGS;
   CFLAGS=; export CFLAGS;  CPPFLAGS=; export CPPFLAGS;

   case "$proglang" in
      Fortran ) AC_LANG([Fortran])     # Fortran (tests assume F90+)
         FCFLAGS="${FCFLAGS_cfg_gen} ${CPPFLAGS_Fortran_cfg_gen}";
         export FCFLAGS;

         FC=${Fortran_compiler};
         AC_PROG_FC

         ftmp=$FCFLAGS;
         AC_FC_SRCEXT([f])
         FCFLAGS=; export FCFLAGS;
         AC_FC_FREEFORM
         FCFLAGS_freeform=$FCFLAGS;
         AC_SUBST([FCFLAGS_freeform])
         FCFLAGS=${ftmp}; export FCFLAGS;

         LDFLAGS="${LDFLAGS_Fortran_cfg_gen}"; export LDFLAGS;
         LIBS="${LIBS_Fortran_cfg_gen}"; export LIBS;

         AC_FC_SRCEXT([f])
         AC_FC_SRCEXT([f95])
         AC_FC_SRCEXT([F])
         AC_FC_SRCEXT([F90])
         AC_FC_SRCEXT([F95])
         AC_FC_SRCEXT([f90])   # Default Fortran extension

         AC_PROG_FC_C_O  # Defines FC_NO_MINUS_C_MINUS_O

         AC_FC_LIBRARY_LDFLAGS

         AC_FC_WRAPPERS

         AC_OPENMP

         # Check compiler for simple Fortran 2003 capabilities:
         AC_LANG_ASSERT([Fortran])
	 AC_MSG_NOTICE([[NOTE: Basic Fortran 2003 capability is necessary for this application]])
         AC_MSG_CHECKING([[whether Fortran compiler has basic Fortran 2003 capabilities]])
         AC_COMPILE_IFELSE(
           [AC_LANG_PROGRAM([],
           [use iso_c_binding, only: C_CHAR, C_NULL_CHAR
            integer :: count
            count = command_argument_count()])],
           [loc_result=yes],[loc_result=no])
         if test "$loc_result" = no; then
            AC_MSG_RESULT([no])
            AC_MSG_FAILURE([[compiler does not seem to have basic Fortran 2003 capabilities -- perhaps add/remove/change a compiler flag?]])
            F2003_ability=no;
         else
            AC_MSG_RESULT([yes])
            F2003_ability=yes;
         fi

         # Save for later:
         LDFLAGS_save="${LDFLAGS}";
         FCFLAGS_save="${FCFLAGS}";
         LIBS_save="${LIBS}";

         for package in $aux_libs_needed
         do
            case $package in
               NetCDF_Fortran )
                  FCFLAGS_to_add=;
                  for tgtdir in ${NetCDF_Fortran_include_dirs}
                  do
                     FCFLAGS_to_add="${FCFLAGS_to_add} -I${tgtdir}";
                  done
                  FCFLAGS="${FCFLAGS_cfg_gen} ${FCFLAGS_to_add}";
                  export FCFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${NetCDF_Fortran_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_Fortran_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${NetCDF_Fortran_libs}";
                  LIBS="${LIBS_Fortran_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  AC_MSG_NOTICE([[NOTE: NetCDF is necessary for this application]])

                  # Check for useable 'netcdf.mod':
                  AC_LANG_ASSERT([Fortran])
                  AC_MSG_CHECKING([[whether NetCDF Fortran 90+ module is present and useable]])

                  AC_COMPILE_IFELSE(
                    [AC_LANG_PROGRAM([],
                    [use netcdf])],
                    [loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_RESULT([no])
                     AC_MSG_FAILURE([[useable netcdf.mod not found -- modify build environment]])
                     failed=yes;
                  else
                     AC_MSG_RESULT([yes])
                  fi

                  if test "$failed" = no; then
                     # Check for useable NetCDF Fortran90+ libraries:
                     AC_LANG_ASSERT([Fortran])
                     AC_MSG_CHECKING([[whether NetCDF Fortran 90+ library is present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([],
                       [[use netcdf, only: nf90_strerror
                         integer :: status
                         status = 1
                         print *,nf90_strerror(status)]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])
                        AC_MSG_FAILURE([[useable libnetcdff not found -- modify build environment]])
                        failed=yes;
                     else
                        AC_MSG_RESULT([yes])
                     fi
                  fi

                  if test "$failed" = no; then
                     FCFLAGS_save="${FCFLAGS_save} ${FCFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;
               HDF4_Fortran )
	          FCFLAGS_to_add=;
                  for tgtdir in ${HDF4_Fortran_include_dirs}
                  do
                     FCFLAGS_to_add="${FCFLAGS_to_add} -I${tgtdir}";
                  done
                  FCFLAGS="${FCFLAGS_cfg_gen} ${FCFLAGS_to_add}";
                  export FCFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${HDF4_Fortran_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_Fortran_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${HDF4_Fortran_libs}";
                  LIBS="${LIBS_Fortran_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  AC_MSG_NOTICE([[NOTE: HDF4 is necessary for this application]])

                  # Check for useable 'hdf.f90':
                  AC_LANG_ASSERT([Fortran])
                  AC_MSG_CHECKING([[whether HDF4 Fortran 90+ include file is present and useable]])
                  AC_COMPILE_IFELSE(
                    [AC_LANG_PROGRAM([],
                    [[include 'hdf.f90'
                      integer :: i
                      i = 7]])],
                    [loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_RESULT([no])
                     AC_MSG_FAILURE([[useable hdf.f90 not found -- modify build environment]])
                     failed=yes;
                  else
                     AC_MSG_RESULT([yes])
                  fi

                  if test "$failed" = no; then
                     # Check for useable HDF4 Fortran 90+ libraries:
                     AC_LANG_ASSERT([Fortran])
                     AC_MSG_CHECKING([[whether HDF4 Fortran 90+ library is present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([],
                       [[include 'hdf.f90'
                         integer :: status
                         status = sfstart('dummy.hdf', DFACC_READ)]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])

                        # Try adding '-ltirpc':
			LDFLAGS_to_add="${LDFLAGS_to_add} -ltirpc";
                        LDFLAGS="${LDFLAGS_Fortran_cfg_gen} ${LDFLAGS_to_add}";
                        export LDFLAGS;
                        AC_LINK_IFELSE(
                          [AC_LANG_PROGRAM([],
                          [[include 'hdf.f90'
                            integer :: status
                            status = sfstart('dummy.hdf', DFACC_READ)]])],
                          [loc_result=yes],[loc_result=no])
                        if test "$loc_result" = no; then
                           AC_MSG_RESULT([no])
                           AC_MSG_FAILURE([[useable libmfhdf not found -- modify build environment]])
                           failed=yes;
                        else
                           AC_MSG_RESULT([yes])
                        fi
                     else
                        AC_MSG_RESULT([yes])
		     fi
                  fi

                  if test "$failed" = no; then
                     FCFLAGS_save="${FCFLAGS_save} ${FCFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;
               HDF5_Fortran )
	          FCFLAGS_to_add=;
                  for tgtdir in ${HDF5_Fortran_include_dirs}
                  do
                     FCFLAGS_to_add="${FCFLAGS_to_add} -I${tgtdir}";
                  done
                  FCFLAGS="${FCFLAGS_cfg_gen} ${FCFLAGS_to_add}";
                  export FCFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${HDF5_Fortran_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_Fortran_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${HDF5_Fortran_libs}";
                  LIBS="${LIBS_Fortran_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  AC_MSG_NOTICE([[NOTE: HDF5 is necessary for this application]])

                  # Check for useable 'hdf5.mod':
                  AC_LANG_ASSERT([Fortran])
                  AC_MSG_CHECKING([[whether HDF5 Fortran 90+ module is present and useable]])

                  AC_COMPILE_IFELSE(
                    [AC_LANG_PROGRAM([],
                    [use hdf5])],
                    [loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_RESULT([no])
                     AC_MSG_FAILURE([[useable hdf5.mod not found -- modify build environment]])
                     failed=yes;
                  else
                     AC_MSG_RESULT([yes])
                  fi

                  if test "$failed" = no; then
                     # Check for useable HDF5 Fortran90+ libraries:
                     AC_LANG_ASSERT([Fortran])
                     AC_MSG_CHECKING([[whether HDF5 Fortran 90+ library is present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([],
                       [[use hdf5, only: h5fopen_f, H5F_ACC_RDONLY_F, HID_T
		         integer(HID_T) :: id_f
                         integer :: status
                         call h5fopen_f('dummy.h5', H5F_ACC_RDONLY_F, id_f, status)]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])
                        AC_MSG_FAILURE([[useable libhdf5_fortran not found -- modify build environment]])
                        failed=yes;
                     else
                        AC_MSG_RESULT([yes])
                     fi
                  fi

                  if test "$failed" = no; then
                     FCFLAGS_save="${FCFLAGS_save} ${FCFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;  
               deGRIB_Fortran )
	          FCFLAGS_to_add=;
                  for tgtdir in ${deGRIB_Fortran_include_dirs}
                  do
                     FCFLAGS_to_add="${FCFLAGS_to_add} -I${tgtdir}";
                  done
                  FCFLAGS="${FCFLAGS_cfg_gen} ${FCFLAGS_to_add}";
                  export FCFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${deGRIB_Fortran_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_Fortran_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${deGRIB_Fortran_libs}";
                  LIBS="${LIBS_Fortran_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  AC_MSG_NOTICE([[NOTE: de-GRIB is an optional capability (can be turned off in build/env_settings/user_change_me.cfg)]]) 

                  # Check for useable 'eccodes.mod':
                  AC_LANG_ASSERT([Fortran])
                  AC_MSG_CHECKING([[whether (de-GRIB) ecCodes Fortran 90+ module is present and useable]])

                  AC_COMPILE_IFELSE(
                    [AC_LANG_PROGRAM([],
                    [use eccodes])],
                    [loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_RESULT([no])
                     AC_MSG_FAILURE([[useable eccodes.mod not found -- modify build environment]])
                     failed=yes;
                  else
                     AC_MSG_RESULT([yes])
                  fi

                  if test "$failed" = no; then
                     # Check for useable (de-GRIB) ecCodes Fortran90+ libraries:
                     AC_LANG_ASSERT([Fortran])
                     AC_MSG_CHECKING([[whether (de-GRIB) ecCodes Fortran 90+ library is present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([],
                       [[use eccodes, only: codes_open_file
                         integer :: status
                         call codes_open_file(status, 'dummy.grib', 'r')]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])
                        AC_MSG_FAILURE([[useable libeccodes_f90 not found -- modify build environment]])
                        failed=yes;
                     else
                        AC_MSG_RESULT([yes])
                     fi
                  fi

                  if test "$failed" = no; then
                     FCFLAGS_save="${FCFLAGS_save} ${FCFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;
	       libHim_Fortran )
	          FCFLAGS_to_add=;
                  for tgtdir in ${libHim_Fortran_include_dirs}
                  do
                     FCFLAGS_to_add="${FCFLAGS_to_add} -I${tgtdir}";
                  done
                  FCFLAGS="${FCFLAGS_cfg_gen} ${FCFLAGS_to_add}";
                  export FCFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${libHim_Fortran_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_Fortran_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${libHim_Fortran_libs}";
                  LIBS="${LIBS_Fortran_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  AC_MSG_NOTICE([[NOTE: using libHimawari is an optional capability (can be turned off in build/env_settings/user_change_me.cfg)]]) 

                  # Check for useable 'eccodes.mod':
                  AC_LANG_ASSERT([Fortran])
                  AC_MSG_CHECKING([[whether libHimawari Fortran 90+ module is present and useable]])

                  AC_COMPILE_IFELSE(
                    [AC_LANG_PROGRAM([],
                    [use eccodes])],
                    [loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_RESULT([no])
                     AC_MSG_FAILURE([[useable eccodes.mod not found -- modify build environment]])
                     failed=yes;
                  else
                     AC_MSG_RESULT([yes])
                  fi

                  if test "$failed" = no; then
                     # Check for useable libHimawari Fortran90+ libraries:
                     AC_LANG_ASSERT([Fortran])
                     AC_MSG_CHECKING([[whether libHimawari Fortran 90+ library is present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([],
                       [[use eccodes, only: codes_open_file
                         integer :: status
                         call codes_open_file(status, 'dummy.grib', 'r')]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])
                        AC_MSG_FAILURE([[useable libeccodes_f90 not found -- modify build environment]])
                        failed=yes;
                     else
                        AC_MSG_RESULT([yes])
                     fi
                  fi

                  if test "$failed" = no; then
                     FCFLAGS_save="${FCFLAGS_save} ${FCFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;
	       RTTOV_Fortran )
	          FCFLAGS_to_add=;
                  for tgtdir in ${RTTOV_Fortran_include_dirs}
                  do
                     FCFLAGS_to_add="${FCFLAGS_to_add} -I${tgtdir}";
                  done
                  FCFLAGS="${FCFLAGS_cfg_gen} ${FCFLAGS_to_add}";
                  export FCFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${RTTOV_Fortran_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_Fortran_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${RTTOV_Fortran_libs}";
                  LIBS="${LIBS_Fortran_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  AC_MSG_NOTICE([[NOTE: using RTTOV is an optional capability (can be turned off in build/env_settings/user_change_me.cfg)]]) 

                  # Check for useable 'eccodes.mod':
                  AC_LANG_ASSERT([Fortran])
                  AC_MSG_CHECKING([[whether RTTOV Fortran 90+ module is present and useable]])

                  AC_COMPILE_IFELSE(
                    [AC_LANG_PROGRAM([],
                    [use eccodes])],
                    [loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_RESULT([no])
                     AC_MSG_FAILURE([[useable eccodes.mod not found -- modify build environment]])
                     failed=yes;
                  else
                     AC_MSG_RESULT([yes])
                  fi

                  if test "$failed" = no; then
                     # Check for useable RTTOV Fortran90+ libraries:
                     AC_LANG_ASSERT([Fortran])
                     AC_MSG_CHECKING([[whether RTTOV Fortran 90+ library is present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([],
                       [[use eccodes, only: codes_open_file
                         integer :: status
                         call codes_open_file(status, 'dummy.grib', 'r')]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])
                        AC_MSG_FAILURE([[useable libeccodes_f90 not found -- modify build environment]])
                        failed=yes;
                     else
                        AC_MSG_RESULT([yes])
                     fi
                  fi

                  if test "$failed" = no; then
                     FCFLAGS_save="${FCFLAGS_save} ${FCFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;
	       CRTM_Fortran )
	          FCFLAGS_to_add=;
                  for tgtdir in ${CRTM_Fortran_include_dirs}
                  do
                     FCFLAGS_to_add="${FCFLAGS_to_add} -I${tgtdir}";
                  done
                  FCFLAGS="${FCFLAGS_cfg_gen} ${FCFLAGS_to_add}";
                  export FCFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${CRTM_Fortran_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_Fortran_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${CRTM_Fortran_libs}";
                  LIBS="${LIBS_Fortran_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  AC_MSG_NOTICE([[NOTE: using CRTM is an optional capability (can be turned off in build/env_settings/user_change_me.cfg)]]) 

                  # Check for useable 'eccodes.mod':
                  AC_LANG_ASSERT([Fortran])
                  AC_MSG_CHECKING([[whether CRTM Fortran 90+ module is present and useable]])

                  AC_COMPILE_IFELSE(
                    [AC_LANG_PROGRAM([],
                    [use eccodes])],
                    [loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_RESULT([no])
                     AC_MSG_FAILURE([[useable eccodes.mod not found -- modify build environment]])
                     failed=yes;
                  else
                     AC_MSG_RESULT([yes])
                  fi

                  if test "$failed" = no; then
                     # Check for useable CRTM Fortran90+ libraries:
                     AC_LANG_ASSERT([Fortran])
                     AC_MSG_CHECKING([[whether CRTM Fortran 90+ library is present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([],
                       [[use eccodes, only: codes_open_file
                         integer :: status
                         call codes_open_file(status, 'dummy.grib', 'r')]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])
                        AC_MSG_FAILURE([[useable libeccodes_f90 not found -- modify build environment]])
                        failed=yes;
                     else
                        AC_MSG_RESULT([yes])
                     fi
                  fi

                  if test "$failed" = no; then
                     FCFLAGS_save="${FCFLAGS_save} ${FCFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;
            esac
         done

         # Separate preprocessor defines and miscellaneous flags from include
         #  directives:

         PREPROC_DEFS_FORTRAN=; PREPROC_INCS_FORTRAN=;
         PREPROC_MFLAGS_FORTRAN=; FLAGS_FORTRAN=;
         for token in $FCFLAGS_save
         do
            if expr X"$token" : '^X\-I' >/dev/null 2>&1; then
               PREPROC_INCS_FORTRAN="$PREPROC_INCS_FORTRAN $token";
            else
               if expr X"$token" : '^X\-D' >/dev/null 2>&1; then
                  PREPROC_DEFS_FORTRAN="$PREPROC_DEFS_FORTRAN $token";
               else
                  FLAGS_FORTRAN="$FLAGS_FORTRAN $token";
               fi
            fi
         done

	 # Finalize and "publish" variables:
         LDFLAGS_FORTRAN="${LDFLAGS_save}";
         LIBS_FORTRAN="${LIBS_save}";
         AC_SUBST([LDFLAGS_FORTRAN])
         AC_SUBST([FLAGS_FORTRAN])
         AC_SUBST([PREPROC_DEFS_FORTRAN])
         AC_SUBST([PREPROC_INCS_FORTRAN])
         AC_SUBST([PREPROC_MFLAGS_FORTRAN])
         AC_SUBST([LIBS_FORTRAN])
           ;;

      C ) AC_LANG([C])        # C (assumes "modern" compiler)
         CFLAGS="${CFLAGS_cfg_gen}"; export CFLAGS;
         CPPFLAGS="${CPPFLAGS_C_cfg_gen}"; export CPPFLAGS;

         CC=${C_compiler};
         AC_PROG_CC                 

         AC_PROG_CPP

         LDFLAGS="${LDFLAGS_C_cfg_gen}"; export LDFLAGS;
         LIBS="${LIBS_C_cfg_gen}"; export LIBS;

         AC_PROG_CC_C_O  # Defines NO_MINUS_C_MINUS

         AC_C_BIGENDIAN([ENDIAN_TYPE=BIG], [ENDIAN_TYPE=LITTLE], [ENDIAN_TYPE=UNIV])
         AC_SUBST([ENDIAN_TYPE])

         AC_C_RESTRICT   # Defines 'restrict'

         AC_C_INLINE   # Defines 'inline'

         AC_C_CHAR_UNSIGNED  # Defines __CHAR_UNSIGNED__

         AC_C_FLEXIBLE_ARRAY_MEMBER  # Defines FLEXIBLE_ARRAY_MEMBER 

         AC_C_VARARRAYS  # Defines HAVE_C_VARARRAYS

         AC_C_TYPEOF  # Defines HAVE_TYPEOF, 'typeof'

         AC_OPENMP

         # Save for later:
         CPPFLAGS_save="${CPPFLAGS}";
         LDFLAGS_save="${LDFLAGS}";
         CFLAGS_save="${CFLAGS}";
         LIBS_save="${LIBS}";

         for package in $aux_libs_needed
         do
            case $package in
               NetCDF_C )
                  CFLAGS_to_add=; CPPFLAGS_to_add=;
                  for tgtdir in ${NetCDF_C_include_dirs}
                  do
                     CFLAGS_to_add="${CFLAGS_to_add} -I${tgtdir}";
                     CPPFLAGS_to_add="${CPPFLAGS_to_add} -I${tgtdir}";
                  done
                  CFLAGS="${CFLAGS_cfg_gen} ${CFLAGS_to_add}";
                  CPPFLAGS="${CPPFLAGS_cfg_gen} ${CPPFLAGS_to_add}";
                  export CFLAGS CPPFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${NetCDF_C_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_C_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${NetCDF_C_libs}";
                  LIBS="${LIBS_C_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  # Check for useable 'netcdf.h':
                  AC_LANG_ASSERT([C])

                  AC_CHECK_HEADERS(netcdf.h,[loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_WARN([[useable netcdf.h not found -- modify build environment]])
                     failed=yes;
                  fi

                  if test "$failed" = no; then
                     # Check for useable NetCDF C libraries:
                     AC_LANG_ASSERT([C])
                     AC_MSG_CHECKING([[whether NetCDF C library is present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([#include <netcdf.h>],
                       [[int ierr;
                         nc_strerror(ierr);]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])
                        AC_MSG_WARN([[useable libnetcdf not found -- modify build environment]])
                        failed=yes;
                     else
                        AC_MSG_RESULT([yes])
                     fi
                  fi

                  if test "$failed" = no; then
                     CPPFLAGS_save="${CPPFLAGS_save} ${CPPFLAGS_to_add}";
                     CFLAGS_save="${CFLAGS_save} ${CFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;
               MPI_C )
                  CFLAGS_to_add=; CPPFLAGS_to_add=;
                  for tgtdir in ${MPI_C_include_dirs}
                  do
                     CFLAGS_to_add="${CFLAGS_to_add} -I${tgtdir}";
                     CPPFLAGS_to_add="${CPPFLAGS_to_add} -I${tgtdir}";
                  done
                  CFLAGS="${CFLAGS_cfg_gen} ${CFLAGS_to_add}";
                  CPPFLAGS="${CPPFLAGS_cfg_gen} ${CPPFLAGS_to_add}";
                  export CFLAGS CPPFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${MPI_C_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_C_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${MPI_C_libs}";
                  LIBS="${LIBS_C_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  # Check for useable 'mpi.h':
                  AC_LANG_ASSERT([C])

                  AC_CHECK_HEADERS(mpi.h,[loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_WARN([[useable mpi.h not found -- modify build environment]])
                     failed=yes;
                  fi

                  if test "$failed" = no; then
                     # Check for useable MPI C libraries:
                     AC_LANG_ASSERT([C])
                     AC_MSG_CHECKING([[whether MPI C library is present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([#include <mpi.h>],
                       [[int ierr;
                         ierr = MPI_Finalize();]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])
                        AC_MSG_WARN([[useable library not found -- modify build environment]])
                        failed=yes;
                     else
                        AC_MSG_RESULT([yes])
                     fi
                  fi

                  if test "$failed" = no; then
                     CPPFLAGS_save="${CPPFLAGS_save} ${CPPFLAGS_to_add}";
                     CFLAGS_save="${CFLAGS_save} ${CFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;
               CFITSIO_C )
                  CFLAGS_to_add=; CPPFLAGS_to_add=;
                  for tgtdir in ${CFITSIO_C_include_dirs}
                  do
                     CFLAGS_to_add="${CFLAGS_to_add} -I${tgtdir}";
                     CPPFLAGS_to_add="${CPPFLAGS_to_add} -I${tgtdir}";
                  done
                  CFLAGS="${CFLAGS_cfg_gen} ${CFLAGS_to_add}";
                  CPPFLAGS="${CPPFLAGS_cfg_gen} ${CPPFLAGS_to_add}";
                  export CFLAGS CPPFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${CFITSIO_C_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_C_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${CFITSIO_C_libs}";
                  LIBS="${LIBS_C_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  # Check for useable 'fitsio.h':
                  AC_LANG_ASSERT([C])

                  AC_CHECK_HEADERS(fitsio.h,[loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_WARN([[useable fitsio.h not found -- modify build environment]])
                     failed=yes;
                  fi

                  if test "$failed" = no; then
                     # Check for useable CFITSIO C libraries:
                     AC_LANG_ASSERT([C])
                     AC_MSG_CHECKING([[whether CFITSIO C library is present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([#include <fitsio.h>],
                       [[int ierr;
                         fits_write_errmark();]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])
                        AC_MSG_WARN([[useable libcfitsio not found -- modify build environment]])
                        failed=yes;
                     else
                        AC_MSG_RESULT([yes])
                     fi
                  fi

                  if test "$failed" = no; then
                     CPPFLAGS_save="${CPPFLAGS_save} ${CPPFLAGS_to_add}";
                     CFLAGS_save="${CFLAGS_save} ${CFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;
               NCL_NCARg_C )
                  CFLAGS_to_add=; CPPFLAGS_to_add=;
                  for tgtdir in ${NCL_NCARg_C_include_dirs}
                  do
                     CFLAGS_to_add="${CFLAGS_to_add} -I${tgtdir}";
                     CPPFLAGS_to_add="${CPPFLAGS_to_add} -I${tgtdir}";
                  done
                  CFLAGS="${CFLAGS_cfg_gen} ${CFLAGS_to_add}";
                  CPPFLAGS="${CPPFLAGS_cfg_gen} ${CPPFLAGS_to_add}";
                  export CFLAGS CPPFLAGS;

                  LDFLAGS_to_add=;
                  for tgtdir in ${NCL_NCARg_C_lib_dirs}
                  do
                     LDFLAGS_to_add="${LDFLAGS_to_add} -L${tgtdir}";
                  done
                  LDFLAGS="${LDFLAGS_C_cfg_gen} ${LDFLAGS_to_add}";
                  export LDFLAGS;

                  LIBS_to_add="${NCL_NCARg_C_libs}";
                  LIBS="${LIBS_C_cfg_gen} ${LIBS_to_add}";
                  export LIBS;

                  failed=no;

                  # Check for useable 'ncargC.h':
                  AC_LANG_ASSERT([C])

                  AC_CHECK_HEADERS(ncarg/ncargC.h,[loc_result=yes],[loc_result=no])
                  if test "$loc_result" = no; then
                     AC_MSG_WARN([[useable ncargC.h not found -- modify build environment]])
                     failed=yes;
                  fi

                  if test "$failed" = no; then
                     # Check for useable NCAR Graphics C libraries:
                     AC_LANG_ASSERT([C])
                     AC_MSG_CHECKING([[whether NCAR Graphics C libraries are present and useable]])
                     AC_LINK_IFELSE(
                       [AC_LANG_PROGRAM([#include <stdio.h>
                         #include <math.h>
                         #include <ncarg/ncargC.h>
                         #include <ncarg/gks.h>],
                       [[gopen_gks ("stdout",0);
                         gclose_gks();]])],
                       [loc_result=yes],[loc_result=no])
                     if test "$loc_result" = no; then
                        AC_MSG_RESULT([no])
                        AC_MSG_WARN([[useable libncarg not found -- modify build environment]])
                        failed=yes;
                     else
                        AC_MSG_RESULT([yes])
                     fi
                  fi

                  if test "$failed" = no; then
                     CPPFLAGS_save="${CPPFLAGS_save} ${CPPFLAGS_to_add}";
                     CFLAGS_save="${CFLAGS_save} ${CFLAGS_to_add}";
                     LDFLAGS_save="${LDFLAGS_save} ${LDFLAGS_to_add}";
                     LIBS_save="${LIBS_save} ${LIBS_to_add}";
                     usable_aux_libs="${usable_aux_libs} ${package}";
                  else
                     failed_aux_libs="${failed_aux_libs} ${package}";
                  fi
               ;;
            esac
         done

         # Separate preprocessor defines and miscellaneous flags from include
         #  directives:

         PREPROC_DEFS_C=; PREPROC_INCS_C=; PREPROC_MFLAGS_C=;
         for token in $CPPFLAGS_save
         do
            if expr X"$token" : '^X\-I' >/dev/null 2>&1; then
               PREPROC_INCS_C="$PREPROC_INCS_C $token";
            else
               if expr X"$token" : '^X\-D' >/dev/null 2>&1; then
                  PREPROC_DEFS_C="$PREPROC_DEFS_C $token";
               else
                  PREPROC_MFLAGS_C="$PREPROC_MFLAGS_C $token";
               fi
            fi
         done

         FLAGS_C=;
         for token in $CFLAGS_save
         do
            if expr X"$token" : '^X\-@<:@^ID@:>@' >/dev/null 2>&1; then
               FLAGS_C="$FLAGS_C $token";
            fi
         done

	 # Finalize and "publish" variables:
         LDFLAGS_C="${LDFLAGS_save}";
         LIBS_C="${LIBS_save}";
         AC_SUBST([LDFLAGS_C])
         AC_SUBST([FLAGS_C])
         AC_SUBST([PREPROC_DEFS_C])
         AC_SUBST([PREPROC_INCS_C])
         AC_SUBST([PREPROC_MFLAGS_C])
         AC_SUBST([LIBS_C])
           ;;
   esac
done

SOME_F2003="${F2003_ability}";
AC_SUBST([SOME_F2003])

FAILED_AUXLIBS="${failed_aux_libs}"; USABLE_AUXLIBS="${usable_aux_libs}";
AC_SUBST([FAILED_AUXLIBS])
AC_SUBST([USABLE_AUXLIBS])

# Output some results to file(s):

AC_OUTPUT
