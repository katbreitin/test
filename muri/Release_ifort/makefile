################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

HOST=$(shell hostname)
ifeq ($(HOST),orchid-submit.ssec.wisc.edu)
LIBS += -lsz -ltirpc
endif
ifeq ($(HOST),vor.ssec.wisc.edu)
LIBS += -lsz -ltirpc
endif
ifeq ($(HOST),saga.ssec.wisc.edu)
LIBS += -lsz -ltirpc
endif
ifeq ($(HOST),solar2.ssec.wisc.edu)
LIBS += -lsz -ltirpc
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: MURI post-build

# Tool invocations
MURI: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: MacOS X Fortran Linker'
	ifort -traceback  -O2 -assume byterecl -profile-functions -profile-loops=all -profile-loops-report=2 -L"${HDF4_PATH}/lib/" -L"${NETCDF_PATH}/lib/" -L"${HDF5_PATH}/lib/" -L"${CX_SDS_IO_PATH}" -o muri_on_pixel $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) MURI
	-@echo ' '

post-build:
	-ar rvc libcx_muri.a src/*.o
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets
