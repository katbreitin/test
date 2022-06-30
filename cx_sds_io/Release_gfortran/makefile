################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include makefile.init

RM := rm -rf

GF_VERSION := $(shell gfortran -dumpversion)
VERSION := 10.
fflags1= -O2
ifeq ($(VERSION),$(firstword $(sort $(VERSION) $(GF_VERSION))))
fflags1 += -fallow-argument-mismatch 
endif

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
all: cx_sds_io




# Tool invocations
cx_sds_io: $(OBJS) $(USER_OBJS)
	
	@echo 'Building target: $@'
	@echo 'Invoking: MacOS X Fortran Linker'
	gfortran $(fflags1) -L"${HDF5_PATH}/lib/" -L"${HDF4_PATH}/lib" -L"${NETCDF_PATH}/lib/" -o "cx_sds_io" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '
	$(MAKE) --no-print-directory post-build

# Other Targets
clean:
	-$(RM) $(EXECUTABLES)$(OBJS)$(C_DEPS) cx_sds_io
	-@echo ' '

post-build:
	-@echo 'Build static library'
	-ar rvc libcx_sds_io.a src/*.o
	-@echo ' '

.PHONY: all clean dependents post-build

-include ../makefile.targets
