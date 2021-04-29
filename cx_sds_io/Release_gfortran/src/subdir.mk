################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../src/cx_h5_read_mod.f90 \
../src/cx_hdf_read_mod.f90 \
../src/cx_hdf_write_mod.f90 \
../src/cx_ncdf_read_mod.f90 \
../src/cx_sds_io_mod.f90 \
../src/cx_sds_test.f90 \
../src/cx_sds_type_definitions_mod.f90 \
../src/cx_sds_write_mod.f90 \
../src/readh5dataset.f90 

OBJS += \
./src/cx_h5_read_mod.o \
./src/cx_hdf_read_mod.o \
./src/cx_hdf_write_mod.o \
./src/cx_ncdf_read_mod.o \
./src/cx_sds_io_mod.o \
./src/cx_sds_test.o \
./src/cx_sds_type_definitions_mod.o \
./src/cx_sds_write_mod.o \
./src/readh5dataset.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -I"${HDF4_PATH}/include" -I"${HDF5_PATH}/include/" -I"${NETCDF_PATH}//include/" -O2  -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/cx_h5_read_mod.o: ../src/cx_h5_read_mod.f90 src/cx_sds_type_definitions_mod.o src/readh5dataset.o

src/cx_hdf_read_mod.o: ../src/cx_hdf_read_mod.f90 src/cx_sds_type_definitions_mod.o

src/cx_hdf_write_mod.o: ../src/cx_hdf_write_mod.f90

src/cx_ncdf_read_mod.o: ../src/cx_ncdf_read_mod.f90 src/cx_sds_type_definitions_mod.o

src/cx_sds_io_mod.o: ../src/cx_sds_io_mod.f90 src/cx_h5_read_mod.o src/cx_hdf_read_mod.o src/cx_ncdf_read_mod.o src/cx_sds_type_definitions_mod.o

src/cx_sds_test.o: ../src/cx_sds_test.f90 src/cx_sds_io_mod.o src/cx_sds_type_definitions_mod.o

src/cx_sds_type_definitions_mod.o: ../src/cx_sds_type_definitions_mod.f90

src/cx_sds_write_mod.o: ../src/cx_sds_write_mod.f90

src/readh5dataset.o: ../src/readh5dataset.f90


