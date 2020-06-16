#
# Build configuration
#

COMP            = gcc
DEBUG           = FALSE
USE_CUDA        = FALSE
USE_KLU         = FALSE
CUDA_ARCH       =

#SS_VERSION=5.7.2
SS_VERSION=5.4.0

ifeq ($(USE_CUDA),TRUE)
  ifeq ($(CUDA_ARCH),)
    $(error Must set CUDA_ARCH.  Run deviceQuery and if CUDA Capability is 6.1, CUDA_ARCH=sm_61. Cori:sm_70, Garuda:sm_75, Godzilla:sm_61)
  endif
endif

lowercase_hcomp = $(shell echo $(COMP) | tr A-Z a-z)
ifeq ($(lowercase_hcomp),$(filter $(lowercase_hcomp),gcc gnu g++ mpicc mpicxx))
  ifeq ($(lowercase_hcomp),$(filter $(lowercase_hcomp),mpi))
    CCOMPILER = $(shell $(COMP) --show | awk '{print $$1}')
  else
    CCOMPILER = gcc
  endif
else
  ifeq ($(lowercase_hcomp),$(filter $(lowercase_hcomp),intel))
    CCOMPILER = icc
  else
    ifeq ($(lowercase_hcomp),$(filter $(lowercase_hcomp),pgi))
      CCOMPILER = pgicc
    else
      $(error Unknown COMP setting)
    endif
  endif
endif

ifeq ($(USE_MPI),TRUE)
  ifeq ($(CCOMPILER),$(filter $(CCOMPILER),mpi))
    HOST_COMPILER = $(shell $(CCOMPILER) --show | awk '{print $$1}')
  else
    HOST_COMPILER = $(CCOMPILER)
  endif
else
  HOST_COMPILER = $(CCOMPILER)
endif

ifeq ($(USE_CUDA),TRUE)
  DEV_COMPILER = nvcc
else
  DEV_COMPILER = $(HOST_COMPILER)
endif

tp_suffix=$(COMP)$(DebugSuffix)$(CUDASuffix)$(KLUSuffix)
INSTALL_PREFIX=$(PWD)/INSTALL/$(tp_suffix)

ifeq ($(USE_KLU),TRUE)
  KLUSuffix   := .KLU
  KLU_DEFINES=-DKLU_ENABLE:BOOL=ON -DKLU_INCLUDE_DIR=${INSTALL_PREFIX}/include -DKLU_LIBRARY_DIR=${INSTALL_PREFIX}/lib 
  SS_BUILD_PREFIX=$(PWD)/BUILD/SUITESPARSE

  SS_DIST_DIR=$(SS_BUILD_PREFIX)/dist
  SS_SOURCE_DIR=$(SS_BUILD_PREFIX)/SuiteSparse-$(SS_VERSION)
  SS_BUILD_DIR=$(SS_BUILD_PREFIX)/build/$(tp_suffix)

  ONE_SS_SRC_FILE=$(SS_SOURCE_DIR)/LICENSE.txt
  ONE_SS_LIB_FILE=$(INSTALL_PREFIX)/libamd.so
  SS_DIST_FILE=v$(SS_VERSION).tar.gz
else
  KLUSuffix=
  KLU_DEFINES=-DKLU_ENABLE:BOOL=OFF
  ONE_SS_SRC_FILE=
  ONE_SS_LIB_FILE=
endif

ifeq ($(DEBUG),FALSE)
  BUILD_TYPE=Release
  COMP_FLAGS="-O3 -DNDEBUG"
  DebugSuffix=
else
  BUILD_TYPE=Debug
  DebugSuffix=.DEBUG
endif

ifeq ($(USE_CUDA),TRUE)
  CUDA_DEFINES=-DCUDA_ENABLE:BOOL=ON -DCUDA_HOST_COMPILER=$(DEV_COMPILER)
  CUDASuffix=.CUDA
else
  CUDA_DEFINES=-DCUDA_ENABLE:BOOL=OFF
  CUDASuffix=
endif

SUNDIALS_BUILD_PREFIX=$(PWD)/BUILD/SUNDIALS
SUNDIALS_SOURCE_DIR=$(SUNDIALS_BUILD_PREFIX)/src
SUNDIALS_BUILD_DIR=$(SUNDIALS_BUILD_PREFIX)/build/$(tp_suffix)
ONE_SUNDIALS_SRC_FILE=$(SUNDIALS_SOURCE_DIR)/config
ONE_SUNDIALS_LIB_FILE=$(INSTALL_PREFIX)/lib/libsundials_cvodes.a

all: $(ONE_SUNDIALS_LIB_FILE)

tp_suffix:
	@echo $(tp_suffix)


$(ONE_SUNDIALS_LIB_FILE): $(ONE_SUNDIALS_SRC_FILE) $(ONE_SS_LIB_FILE)
	@mkdir -p $(SUNDIALS_BUILD_DIR)
	cd $(SUNDIALS_BUILD_DIR); cmake -DCUDA_ARCH:STRING=$(CUDA_ARCH) -DCMAKE_INSTALL_PREFIX:STRING=$(INSTALL_PREFIX) -DCMAKE_INSTALL_LIBDIR:STRING=$(INSTALL_PREFIX)/lib -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DCMAKE_C_COMPILER:STRING=$(HOST_COMPILER) $(CUDA_HOST_COMP_DEFINE) -DEXAMPLES_INSTALL_PATH:STRING=$(INSTALL_PREFIX)/examples -DCMAKE_BUILD_TYPE:STRING=$(BUILD_TYPE) -DCMAKE_C_FLAGS_RELEASE:STRING=$(COMP_FLAGS) $(CUDA_DEFINES) -DMPI_ENABLE:BOOL=OFF -DOPENMP_ENABLE:BOOL=OFF $(KLU_DEFINES) -DEXAMPLES_ENABLE_C:BOOL=OFF -DSUNDIALS_INDEX_SIZE:INT=32 $(SUNDIALS_SOURCE_DIR); make; make install

$(ONE_SUNDIALS_SRC_FILE):
	@mkdir -p $(SUNDIALS_SOURCE_DIR)
	@git clone https://github.com/LLNL/sundials $(SUNDIALS_SOURCE_DIR)

$(ONE_SS_LIB_FILE): $(ONE_SS_SRC_FILE)
	@cd $(SS_SOURCE_DIR)/SuiteSparse_config; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/AMD; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/COLAMD; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/BTF; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/metis-5.1.0; make config BLAS=-lblas prefix=$(INSTALL_PREFIX) shared=1; make install
	@cd $(SS_SOURCE_DIR)/CCOLAMD; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/CAMD; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/KLU; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)

$(ONE_SS_SRC_FILE): $(SS_DIST_DIR)/$(SS_DIST_FILE)
	@cd $(SS_DIST_DIR)/..; tar zxf $(SS_DIST_DIR)/$(SS_DIST_FILE); find $(SS_SOURCE_DIR) -type f -exec touch {} +

$(SS_DIST_DIR)/$(SS_DIST_FILE):
	@mkdir -p $(SS_DIST_DIR)
	@cd $(SS_DIST_DIR); wget https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/$(SS_DIST_FILE)

ifeq ($(USE_KLU),TRUE)
  clean_SS:
	@echo rm -rf $(SS_BUILD_PREFIX)
else
  clean_SS:
endif

clean: clean_SS
	rm -rf $(INSTALL_PREFIX) $(SUNDIALS_BUILD_PREFIX)

realclean:
	rm -rf $(PWD)/INSTALL $(PWD)/BUILD
