

# default for compiler architecture macro is gfortran
# any value given on the command line gets higher priority
ARCH=gfortran

# define local flags in MAKEFILE_LOCAL
#

#MAKEFILE_LOCAL=$(TOPDIR)/build/Makefile.local
MAKEFILE_LOCAL=$(TOPDIR)/build/Makefile.$(HOST)
include $(MAKEFILE_LOCAL)


include $(TOPDIR)/build/arch/$(ARCH)

export CMASK=../cloud_mask/
export BASELINEMASK=../baseline_cloud_mask/
export CTYPE=../cloud_type/
export EXTRA=../extra_src/
export ACHA=../acha/
export CCL=../ccl/
export ASOS=../asos/
export PFAAST=../pfaast/
export SASRAB=../sasrab/
export BASE=../cloud_base/
export BASELINEHEIGHT=../baseline_cloud_height/
dncomplibs=-L$(DNCOMP) -I$(DNCOMP) -L$(CX_DATA_IO) -I$(CX_DATA_IO) 
dncomplinks= -lcx_dncomp -lcx_sds_io
muri_libs=-L$(MURI) -I$(MURI)  -lcx_muri
murilinks=-lcx_muri
ifdef RTTOV_PATH
    rttov_libs=-L$(RTTOV_DIR) -I$(RTTOV_DIR)  -lcx_rttov
    rttovlinks=-lcx_rttov    
endif

ifdef LIBHIM_PATH
    libhim_libs=-I$(LIBHIM_PATH)/include -I$(LIBHIM_PATH)/src
    libhim=-lstdc++ $(LIBHIM_PATH)/src/libHimawari.a 
endif
