################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../src/M_kracken.f90 \
../src/dcomp_array_loop_sub.f90 \
../src/dcomp_forward_mod.f90 \
../src/dcomp_lut_mod.f90 \
../src/dcomp_lut_sds_mod.f90 \
../src/dcomp_math_tools_mod.f90 \
../src/dcomp_one_pixel_run_program.f90 \
../src/dcomp_pixel_input__define_mod.f90 \
../src/dcomp_retrieval_mod.f90 \
../src/dcomp_science_tools_mod.f90 \
../src/dncomp_interface_def_mod.f90 \
../src/dncomp_precip_mod.f90 \
../src/dncomp_trans_atmos_mod.f90 \
../src/file_tools.f90 \
../src/nlcomp_array_loop_sub.f90 \
../src/nlcomp_forward_mod.f90 \
../src/nlcomp_retrieval_mod.f90 \
../src/sensorname_from_wmoid.f90 

OBJS += \
./src/M_kracken.o \
./src/dcomp_array_loop_sub.o \
./src/dcomp_forward_mod.o \
./src/dcomp_lut_mod.o \
./src/dcomp_lut_sds_mod.o \
./src/dcomp_math_tools_mod.o \
./src/dcomp_one_pixel_run_program.o \
./src/dcomp_pixel_input__define_mod.o \
./src/dcomp_retrieval_mod.o \
./src/dcomp_science_tools_mod.o \
./src/dncomp_interface_def_mod.o \
./src/dncomp_precip_mod.o \
./src/dncomp_trans_atmos_mod.o \
./src/file_tools.o \
./src/nlcomp_array_loop_sub.o \
./src/nlcomp_forward_mod.o \
./src/nlcomp_retrieval_mod.o \
./src/sensorname_from_wmoid.o 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	gfortran -funderscoring -I"${HDF4_PATH}/include" -I"${HDF5_PATH}/include/" -I"${NETCDF_PATH}/include/" -I"${CX_SDS_IO_PATH}" -O2  -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/M_kracken.o: ../src/M_kracken.f90

src/dcomp_array_loop_sub.o: ../src/dcomp_array_loop_sub.f90 src/dcomp_retrieval_mod.o src/dncomp_interface_def_mod.o src/dncomp_precip_mod.o src/dncomp_trans_atmos_mod.o

src/dcomp_forward_mod.o: ../src/dcomp_forward_mod.f90 src/dcomp_lut_mod.o src/dcomp_science_tools_mod.o

src/dcomp_lut_mod.o: ../src/dcomp_lut_mod.f90 src/dcomp_lut_sds_mod.o src/dcomp_math_tools_mod.o src/file_tools.o

src/dcomp_lut_sds_mod.o: ../src/dcomp_lut_sds_mod.f90

src/dcomp_math_tools_mod.o: ../src/dcomp_math_tools_mod.f90

src/dcomp_one_pixel_run_program.o: ../src/dcomp_one_pixel_run_program.f90 src/M_kracken.o src/dcomp_lut_mod.o src/dcomp_math_tools_mod.o src/dcomp_retrieval_mod.o src/nlcomp_retrieval_mod.o

src/dcomp_pixel_input__define_mod.o: ../src/dcomp_pixel_input__define_mod.f90

src/dcomp_retrieval_mod.o: ../src/dcomp_retrieval_mod.f90 src/dcomp_forward_mod.o src/dcomp_math_tools_mod.o

src/dcomp_science_tools_mod.o: ../src/dcomp_science_tools_mod.f90

src/dncomp_interface_def_mod.o: ../src/dncomp_interface_def_mod.f90

src/dncomp_precip_mod.o: ../src/dncomp_precip_mod.f90

src/dncomp_trans_atmos_mod.o: ../src/dncomp_trans_atmos_mod.f90

src/file_tools.o: ../src/file_tools.f90

src/nlcomp_array_loop_sub.o: ../src/nlcomp_array_loop_sub.f90 src/dncomp_interface_def_mod.o src/dncomp_precip_mod.o src/nlcomp_retrieval_mod.o

src/nlcomp_forward_mod.o: ../src/nlcomp_forward_mod.f90 src/dcomp_lut_mod.o src/dcomp_science_tools_mod.o

src/nlcomp_retrieval_mod.o: ../src/nlcomp_retrieval_mod.f90 src/dcomp_math_tools_mod.o src/nlcomp_forward_mod.o

src/sensorname_from_wmoid.o: ../src/sensorname_from_wmoid.f90


