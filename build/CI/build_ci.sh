#!/bin/bash

export PATH=/usr/local/gcc-8.3.0/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/gcc-8.3.0/lib64

export SSEC_HDF4_INC=/root/hdf4/include
export SSEC_HDF4_LIB=/root/hdf4/lib
export SSEC_HDF5_INC=/root/hdf5/include
export SSEC_HDF5_LIB=/root/hdf5/lib
export SSEC_NETCDF4_LIB=/root/netcdf/lib
export SSEC_NETCDF4_INC=/root/netcdf/include
#export LIBHIM_PATH=/host/out/himawari
#export RTTOV_PATH=/root/rttov/

cd build
cp CI/user_change_me.cfg env_settings
./admin clean ALL
./admin build
