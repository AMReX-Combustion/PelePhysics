COMP            ?= gcc
DEBUG           ?= FALSE
USE_CUDA        ?= FALSE
USE_KLU_PP      ?= FALSE
CUDA_ARCH       ?= 70

SS_VERSION=5.4.0

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

tp_suffix=$(COMP)$(DebugSuffix)$(CUDASuffix)$(KLUSuffix)
INSTALL_PREFIX=$(PWD)/INSTALL/$(tp_suffix)

ifeq ($(USE_KLU_PP),TRUE)
  KLUSuffix   := .KLU
  KLU_DEFINES=-DENABLE_KLU:BOOL=ON -DKLU_INCLUDE_DIR=${INSTALL_PREFIX}/include -DKLU_LIBRARY_DIR=${INSTALL_PREFIX}/lib 
  SS_BUILD_PREFIX=$(PWD)/BUILD/SUITESPARSE

  SS_DIST_DIR=$(SS_BUILD_PREFIX)/dist
  SS_SOURCE_DIR=$(SS_BUILD_PREFIX)/SuiteSparse-$(SS_VERSION)
  SS_BUILD_DIR=$(SS_BUILD_PREFIX)/build/$(tp_suffix)

  ONE_SS_SRC_FILE=$(SS_SOURCE_DIR)/LICENSE.txt
  ONE_SS_LIB_FILE=$(INSTALL_PREFIX)/libamd.so
  SS_DIST_FILE=v$(SS_VERSION).tar.gz
else
  KLUSuffix=
  KLU_DEFINES=-DENABLE_KLU:BOOL=OFF
  ONE_SS_SRC_FILE=
  ONE_SS_LIB_FILE=
endif

ifeq ($(DEBUG),FALSE)
  BUILD_TYPE=Release
  DebugSuffix=
else
  BUILD_TYPE=Debug
  DebugSuffix=.DEBUG
endif

ifeq ($(USE_CUDA),TRUE)
  CUDA_DEFINES=-DENABLE_CUDA:BOOL=ON -DCMAKE_CUDA_ARCHITECTURES:STRING=$(CUDA_ARCH)
  CXX=$(NVCC_CCBIN)
  CC=$(NVCC_CCBIN)
  CUDASuffix=.CUDA
else
  CUDA_DEFINES=-DENABLE_CUDA:BOOL=OFF
  CUDASuffix=
endif
ifeq ($(USE_HIP),TRUE)
  HIP_DEFINES=-DENABLE_HIP:BOOL=ON
  HIPSuffix=.HIP
else
  HIP_DEFINES=-DENABLE_HIP:BOOL=OFF
  HIPSuffix=
endif

SUNDIALS_BUILD_PREFIX=$(PWD)/BUILD/SUNDIALS
SUNDIALS_SOURCE_DIR=$(SUNDIALS_BUILD_PREFIX)/src
SUNDIALS_BUILD_DIR=$(SUNDIALS_BUILD_PREFIX)/build/$(tp_suffix)
ONE_SUNDIALS_SRC_FILE=$(SUNDIALS_SOURCE_DIR)/include/sundials/sundials_config.in
ONE_SUNDIALS_LIB_FILE=$(INSTALL_PREFIX)/lib/libsundials_cvode.a

all: $(ONE_SUNDIALS_LIB_FILE)

tp_suffix:
	@echo $(tp_suffix)


$(ONE_SUNDIALS_LIB_FILE): $(ONE_SUNDIALS_SRC_FILE) $(ONE_SS_LIB_FILE)
	@mkdir -p $(SUNDIALS_BUILD_DIR)
	cd $(SUNDIALS_BUILD_DIR); cmake -DCMAKE_BUILD_WITH_INSTALL_RPATH:BOOL=ON -DCMAKE_INSTALL_PREFIX:STRING=$(INSTALL_PREFIX) \
			   		-DCMAKE_INSTALL_LIBDIR:STRING=$(INSTALL_PREFIX)/lib -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
					-DCMAKE_C_COMPILER:STRING=$(CC) -DCMAKE_CXX_COMPILER:STRING=$(CXX) \
                                        -DEXAMPLES_INSTALL_PATH:STRING=$(INSTALL_PREFIX)/examples \
                                        -DCMAKE_BUILD_TYPE:STRING=$(BUILD_TYPE) $(CUDA_DEFINES) $(HIP_DEFINES) \
                                        -DENABLE_MPI:BOOL=OFF -DENABLE_OPENMP:BOOL=OFF $(KLU_DEFINES) \
					-DEXAMPLES_ENABLE_C:BOOL=OFF -DEXAMPLES_ENABLE_CXX:BOOL=OFF -DEXAMPLES_ENABLE_CUDA:BOOL=OFF \
					-DBUILD_CVODES:BOOL=OFF -DBUILD_IDA:BOOL=OFF -DBUILD_IDAS:BOOL=OFF -DBUILD_KINSOL:BOOL=OFF \
					-DSUNDIALS_INDEX_SIZE:INT=32 -DSUNDIALS_RAJA_BACKENDS:STRING=CUDA $(SUNDIALS_SOURCE_DIR); \
					make; make install

$(ONE_SUNDIALS_SRC_FILE):
	@mkdir -p $(SUNDIALS_SOURCE_DIR)
	@[ "$(shell ls -A $(SUNDIALS_SOURCE_DIR))" ] || git clone https://github.com/LLNL/sundials $(SUNDIALS_SOURCE_DIR)

$(ONE_SS_LIB_FILE): $(ONE_SS_SRC_FILE)
	@cd $(SS_SOURCE_DIR)/SuiteSparse_config; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/AMD; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/COLAMD; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/BTF; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/metis-5.1.0; make config BLAS=-lblas prefix=$(INSTALL_PREFIX) shared=1; make install
	@cd $(SS_SOURCE_DIR)/CCOLAMD; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/CAMD; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	@cd $(SS_SOURCE_DIR)/KLU; make install BLAS=-lblas INSTALL=$(INSTALL_PREFIX)
	# For some reason, the static libraries are built but not installed. Add
	# a workaround here for now, but this should be fixed in SuiteSparse.
	@install -m 644 $(SS_SOURCE_DIR)/SuiteSparse_config/libsuitesparseconfig.a $(INSTALL_PREFIX)/lib
	@install -m 644 $(SS_SOURCE_DIR)/*/Lib/*.a $(INSTALL_PREFIX)/lib

$(ONE_SS_SRC_FILE): $(SS_DIST_DIR)/$(SS_DIST_FILE)
	@cd $(SS_DIST_DIR)/..; tar zxf $(SS_DIST_DIR)/$(SS_DIST_FILE); find $(SS_SOURCE_DIR) -type f -exec touch {} +

$(SS_DIST_DIR)/$(SS_DIST_FILE):
	@mkdir -p $(SS_DIST_DIR)
	@cd $(SS_DIST_DIR); wget https://github.com/DrTimothyAldenDavis/SuiteSparse/archive/$(SS_DIST_FILE)

ifeq ($(USE_KLU_PP),TRUE)
  clean_SS:
	@echo rm -rf $(SS_BUILD_PREFIX)
else
  clean_SS:
endif

clean: clean_SS
	rm -rf $(INSTALL_PREFIX) $(SUNDIALS_BUILD_PREFIX)

realclean:
	rm -rf $(PWD)/INSTALL $(PWD)/BUILD
#-----------------------------------------------------------------------------
# for debugging.  To see the value of a Makefile variable,
# e.g. libraries, simply do "make print-libraries".  This will
# print out the value.
print-%:
	@echo $* is '$($*)'
	@echo '    origin = $(origin $*)'
	@echo '     value = $(value  $*)'
