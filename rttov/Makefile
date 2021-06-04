RTTOVDIR=/Users/awalther/SOFTWARE/RTM/RTTOV123/ 
fc = gfortran -I$(RTTOVDIR)/mod -I$(RTTOVDIR)/include -fPIC -O3 -fopenmp -ffree-line-length-none -cpp 
pfaast_dir = ../





.SUFFIXES: .f90 .F .o

RTTOVDIR=/Users/awalther/SOFTWARE/RTM/RTTOV123/

FC=gfortran
FFLAGS= -I$(RTTOVDIR)/mod -I$(RTTOVDIR)/include -fPIC -O3 -fopenmp -ffree-line-length-none -cpp -DRTTOVPATH=\"$(RTTOVDIR)\"
LDFLAGS_ARCH=-fopenmp

LDFLAGS_HDF5=-L$(SSEC_HDF5_LIB)/../lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz

RTTOV_VERSION=12
EXT_BIN=.exe
EXT_OBJ=.o
EXT_LIB=.a
OBJDIR=$(RTTOVDIR)/obj
BINDIR=$(RTTOVDIR)/bin


#### Do not edit beyond this line ####

#include ../../build/Makefile.local

LIBDIR=$(RTTOVDIR)/lib

ifdef LDFLAGS_HDF5
	LIBS=\
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_brdf_atlas$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_emis_atlas$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_mw_scatt$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_other$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_hdf$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_parallel$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

	LLIBS=-L$(LIBDIR) \
	  -lrttov$(RTTOV_VERSION)_brdf_atlas \
	  -lrttov$(RTTOV_VERSION)_emis_atlas \
	  -lrttov$(RTTOV_VERSION)_mw_scatt \
	  -lrttov$(RTTOV_VERSION)_other \
	  -lrttov$(RTTOV_VERSION)_coef_io \
	  -lrttov$(RTTOV_VERSION)_hdf \
	  -lrttov$(RTTOV_VERSION)_parallel \
	  -lrttov$(RTTOV_VERSION)_main
else
	LIBS=\
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_brdf_atlas$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_emis_atlas$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_mw_scatt$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_other$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_coef_io$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_parallel$(EXT_LIB) \
	  $(LIBDIR)/librttov$(RTTOV_VERSION)_main$(EXT_LIB)

	LLIBS=-L$(LIBDIR) \
	  -lrttov$(RTTOV_VERSION)_brdf_atlas \
	  -lrttov$(RTTOV_VERSION)_emis_atlas \
	  -lrttov$(RTTOV_VERSION)_mw_scatt \
	  -lrttov$(RTTOV_VERSION)_other \
	  -lrttov$(RTTOV_VERSION)_coef_io \
	  -lrttov$(RTTOV_VERSION)_parallel \
	  -lrttov$(RTTOV_VERSION)_main
endif
 
#Compilation rule for f90 files
.f90.o:
	$(fc) -c  $<
   
OBJS = cx_rttov_constants_mod.o \
cx_rttov_bridge_mod.o \
cx_rttov_test.o \
cx_rttov_mapping_mod.o

.PHONY: cx_rttov clean
  

  
all: $(OBJS)
	@echo 'Building target: $@'
	$(fc) -o "test_rttov" $(OBJS)    \
	$(LLIBS)  \
	$(LDFLAGS_NETCDF) $(LDFLAGS_HDF5) $(LDFLAGS_ARCH) $(LDFLAGS_LAPACK)

cx_rttov_test$(EXT_OBJ)::cx_rttov_test.f90 
		$(FC) $(FFLAGS) -c cx_rttov_test.f90 -o cx_rttov_test$(EXT_OBJ)

cx_rttov_test$(EXT_BIN)::cx_rttov_test$(EXT_OBJ) $(LIBS)
	$(fc) -o cx_rttov_test$(EXT_BIN) \
  	cx_rttov_test$(EXT_OBJ) \
	$(LLIBS)  \
	$(LDFLAGS_NETCDF) $(LDFLAGS_HDF5) $(LDFLAGS_ARCH) $(LDFLAGS_LAPACK)	
  
  
cx_rttov_constants_mod$(EXT_OBJ)::cx_rttov_constants_mod.f90 
	$(fc) -c cx_rttov_constants_mod.f90      


cx_rttov_bridge_mod$(EXT_OBJ):: cx_rttov_bridge_mod.f90  cx_rttov_constants_mod.o cx_rttov_mapping_mod.o
	$(FC) $(FFLAGS) -c cx_rttov_bridge_mod.f90  \
	$(LLIBS)  \
	$(LDFLAGS_NETCDF) $(LDFLAGS_HDF5) $(LDFLAGS_ARCH) $(LDFLAGS_LAPACK)	

cx_rttov_mapping_mod$(EXT_OBJ)::cx_rttov_mapping_mod.f90 
		$(FC) $(FFLAGS) -c cx_rttov_mapping_mod.f90


clean:
	rm -f *.o *.mod
	rm -f test_rttov
