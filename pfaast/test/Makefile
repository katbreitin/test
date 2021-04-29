RTTOVDIR=/Users/awalther/SOFTWARE/RTM/RTTOV123/
fc = gfortran -g -O0
fc = gfortran -I$(RTTOVDIR)/mod -I$(RTTOVDIR)/include -fPIC -O3 -fopenmp -ffree-line-length-none -cpp 
main_dir=../../main_src/
pfaast_dir = ../





.SUFFIXES: .f90 .F .o

RTTOVDIR=/Users/awalther/SOFTWARE/RTM/RTTOV123/

FC=gfortran
FFLAGS= -I$(RTTOVDIR)/mod -I$(RTTOVDIR)/include -fPIC -O3 -fopenmp -ffree-line-length-none -cpp
LDFLAGS_ARCH=-fopenmp


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
   
OBJS = cx_pfaast_coef_mod.o \
cx_pfaast_test.o \
cx_string_tools_mod.o 

OBJS_RTM = cx_pfaast_coef_mod.o \
cx_pfaast_rtm_test.o \
cx_string_tools_mod.o \
cx_pfaast_mod.o \
cx_pfaast_tools_mod.o \
cx_pfaast_constants_mod.o

OBJS_RTTOV = cx_pfaast_coef_mod.o \
cx_pfaast_rttov_test.o \
cx_string_tools_mod.o \
cx_pfaast_mod.o \
cx_pfaast_tools_mod.o \
cx_pfaast_constants_mod.o \
cx_rttov_bridge_mod.o

.PHONY: test clean

test: $(OBJS)
	@echo 'Building target: $@'
	$(fc) -o "test_pfaast" $(OBJS) -I../main_src/

test_rtm: $(OBJS_RTM)
	@echo 'Building target: $@'
	$(fc) -o "test_rtm" $(OBJS_RTM) -I../main_src/
  
test_rttov: $(OBJS_RTTOV)
	@echo 'Building target: $@'
	$(FC) -o "test_rttov" $(OBJS_RTTOV) -I../main_src/   \
	$(LLIBS)  \
	$(LDFLAGS_NETCDF) $(LDFLAGS_HDF5) $(LDFLAGS_ARCH) $(LDFLAGS_LAPACK)
  

cx_string_tools_mod$(EXT_OBJ):: $(main_dir)/cx_string_tools_mod.f90  
	$(fc) -c $(main_dir)/cx_string_tools_mod.f90
cx_pfaast_coef_mod$(EXT_OBJ):: $(pfaast_dir)/cx_pfaast_coef_mod.f90 cx_string_tools_mod.o
	$(fc) -c $(pfaast_dir)/cx_pfaast_coef_mod.f90 cx_string_tools_mod.o
cx_pfaast_test$(EXT_OBJ)::cx_pfaast_test.f90 cx_pfaast_coef_mod.o
cx_pfaast_mod$(EXT_OBJ)::$(pfaast_dir)/cx_pfaast_mod.f90 cx_pfaast_tools_mod.o cx_pfaast_constants_mod.o
	$(fc) -c $(pfaast_dir)/cx_pfaast_mod.f90 cx_pfaast_tools_mod.o cx_pfaast_constants_mod.o
cx_pfaast_tools_mod$(EXT_OBJ)::$(pfaast_dir)/cx_pfaast_tools_mod.f90 cx_pfaast_constants_mod.o
	$(fc) -c $(pfaast_dir)/cx_pfaast_tools_mod.f90 cx_pfaast_constants_mod.o
cx_pfaast_constants_mod$(EXT_OBJ)::$(pfaast_dir)/cx_pfaast_constants_mod.f90 
	$(fc) -c $(pfaast_dir)/cx_pfaast_constants_mod.f90      
cx_pfaast_rtm_test$(EXT_OBJ)::cx_pfaast_rtm_test.f90 cx_pfaast_mod.o cx_pfaast_constants_mod.o


cx_pfaast_rttov_test$(EXT_OBJ): cx_pfaast_rttov_test.f90 cx_pfaast_mod.o cx_pfaast_constants_mod.o cx_rttov_bridge_mod.o
	$(FC) $(FFLAGS) -c cx_pfaast_rttov_test.f90 -o cx_pfaast_rttov_test$(EXT_OBJ)

cx_pfaast_rttov_test$(EXT_BIN): cx_pfaast_rttov_test$(EXT_OBJ) $(LIBS)
	$(FC) -o cx_pfaast_rttov_test$(EXT_BIN) \
	cx_pfaast_rttov_test$(EXT_OBJ) \
	$(LLIBS)  \
	$(LDFLAGS_NETCDF) $(LDFLAGS_HDF5) $(LDFLAGS_ARCH) $(LDFLAGS_LAPACK)	

cx_rttov_bridge_mod$(EXT_OBJ): cx_rttov_bridge_mod.f90  cx_pfaast_constants_mod.o
	$(FC) $(FFLAGS) -c cx_rttov_bridge_mod.f90 -o cx_rttov_bridge_mod$(EXT_OBJ) \
	$(LLIBS)  \
	$(LDFLAGS_NETCDF) $(LDFLAGS_HDF5) $(LDFLAGS_ARCH) $(LDFLAGS_LAPACK)	




clean:
	rm -f *.o *.mod
	rm test_pfaast
