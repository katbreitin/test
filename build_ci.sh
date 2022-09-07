#!/bin/bash
export PATH=/usr/local/gcc-8.3.0/bin/:$PATH
export LD_LIBRARY_PATH=/usr/local/gcc-8.3.0/lib64/
export SSEC_HDF4_LIB=/root/hdf4/lib/
export SSEC_HDF5_LIB=/root/hdf5/lib/
export SSEC_NETCDF4_DIR=/root/netcdf
export SSEC_NETCDF4_LIB=/root/netcdf/lib
export SSEC_NETCDF4_INC=/root/netcdf/include
#export LIBHIM_PATH=/host/out/himawari
#export RTTOV_PATH=/root/rttov/
export FC=gfortran
export FFLAGS=-g
export fflags=-g

export ARCH=gfortran
#./configure -LIBHIM_GCC=$LIBHIM_PATH --rttov_path=/host/out/rttov
./configure
cd main_src
export FFLAGS=-g
export fflags=-g
make all_plus
#bash /host/relink.sh
#cp /root/clavrx_trunk/clavrx_bin/clavrxorb /host/out/
