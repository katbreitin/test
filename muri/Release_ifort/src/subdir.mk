################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../src/aw_lib_array.f90 \
../src/lib_array.f90 \
../src/muri_array_loop_sub.f90 \
../src/muri_definitions_mod.f90 \
../src/muri_forward_mod.f90 \
../src/muri_interface_mod.f90 \
../src/muri_lut_mod.f90 \
../src/muri_one_pixel_run.f90 \
../src/muri_atmospheric_corr_mod.f90 \
../src/muri_land_lut_mod.f90 \
../src/muri_land_retrieval_mod.f90 \
../src/muri_retrieval_mod.f90 

OBJS += \
./src/aw_lib_array.o \
./src/lib_array.o \
./src/muri_array_loop_sub.o \
./src/muri_definitions_mod.o \
./src/muri_forward_mod.o \
./src/muri_interface_mod.o \
./src/muri_lut_mod.o \
./src/muri_one_pixel_run.o \
./src/muri_atmospheric_corr_mod.o \
./src/muri_land_retrieval_mod.o \
./src/muri_land_lut_mod.o \
./src/muri_retrieval_mod.o 


# Each subdirectory must supply rules for building srcs it contributes
src/%.o: ../src/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	ifort -I"${HDF4_PATH}/include" -I"${HDF5_PATH}/include/" -I"${NETCDF_PATH}/include/" -I"${CX_SDS_IO_PATH}"  -O2  -c -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/aw_lib_array.o: ../src/aw_lib_array.f90 src/lib_array.o

src/muri_land_retrieval_mod.o: ../src/muri_land_retrieval_mod.f90 src/muri_land_lut_mod.o

src/lib_array.o: ../src/lib_array.f90

src/muri_land_lut_mod.o: ../src/muri_land_lut_mod.f90

src/muri_array_loop_sub.o: ../src/muri_array_loop_sub.f90 src/muri_definitions_mod.o src/muri_interface_mod.o src/muri_retrieval_mod.o src/muri_atmospheric_corr_mod.o src/muri_land_retrieval_mod.o

src/muri_definitions_mod.o: ../src/muri_definitions_mod.f90

src/muri_forward_mod.o: ../src/muri_forward_mod.f90 src/muri_lut_mod.o

src/muri_interface_mod.o: ../src/muri_interface_mod.f90

src/muri_lut_mod.o: ../src/muri_lut_mod.f90 src/aw_lib_array.o

src/muri_one_pixel_run.o: ../src/muri_one_pixel_run.f90 src/muri_definitions_mod.o src/muri_retrieval_mod.o src/muri_atmospheric_corr_mod.o src/muri_land_retrieval_mod.o

src/muri_retrieval_mod.o: ../src/muri_retrieval_mod.f90 src/lib_array.o src/muri_definitions_mod.o src/muri_lut_mod.o

src/muri_atmospheric_corr_mod.o: ../src/muri_atmospheric_corr_mod.f90


