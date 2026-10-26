MODEL             = model600
FC                = mpiifx
CC                = mpiicx
CXX               = mpiicpx
 
FFLAGS += -mcmodel=medium -O2
FFLAGS += -g -traceback -check bounds #debug Chizhou
 
FFLAGS += -march=skylake-avx512 -axCORE-AVX512 # for AMD
# FFLAGS += -xHost                             # for Intel
# FFLAGS += -march=sapphirerapids
# FFLAGS += -double-size 64 -integer-size 32 #-DFUNNELED
# FFLAGS += -fpp -mcmodel=medium -axCORE-AVX512 -align array64byte -qopt-zmm-usage=high
 
EXTRA_FLAGS += -lstdc++
COMPILER_FAMILY   = intel
# Never use STRUMPACK and PASTIC at the same time 
USE_PASTIX        = 0
USE_PASTIX_MURGE  = 0
USE_PASTIX6       = 0
USE_STRUMPACK     = 1
USE_MUMPS         = 0
USE_BLOCK         = 1
USE_HDF5          = 1
USE_FFTW          = 1
USE_MKL           = 0
USE_DIRECT_CONSTRUCTION=0
USE_BICGSTAB      = 0
USE_INTSIZE64     = 0
DEBUG             = 0
LIBS += -liomp5 -pthread -ldl -lm
 
# TODO: MUMPS for USE_MUMPS needs to be compiled
# TODO: 64bits libraries for USE_INTSIZE64 need to be compiled
# LIB64DIR =
 
LIBDIR = /pitagora_work/FUSIO_BA7DISRE/JOREK_LIBRARIES/libraries_recompile_intel
#LIBDIR = /pitagora_work/FUPA1_MHD/JOREK_LIBRARIES/libraries_recompile_amd
# LIBDIR = /pitagora_work/FUPA1_MHD/JOREK_LIBRARIES/libraries_recompile_intel
 
ifeq (1, $(USE_PASTIX))
 
 MKLLIB += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
 MKLINC = -I${MKLROOT}/include
 PASTIX_DIR = $(LIBDIR)/pastix_5.2.3/install
 SCOTCH_DIR = $(LIBDIR)/scotch_5.1.12
 ifeq (1, $(USE_INTSIZE64))
   PASTIX_DIR = $(LIBDIR)/pastix_5.2.3_i64/install
   SCOTCH_DIR = $(LIBDIR)/scotch_5.1.12_i64
 endif
 
 LIB_PASTIX = -L$(PASTIX_DIR) -lpastix -L$(SCOTCH_DIR)/lib -lscotch -lscotcherr
 INC_PASTIX = -I$(PASTIX_DIR)/include
 
 INC_PASTIX += $(MKLINC)
 LIB_PASTIX += $(MKLLIB)
 
endif
 
ifeq (1, $(USE_MUMPS))
# --- MUMPS
 MKLLIB += -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKLROOT}/include
 MUMPS_DIR = $(LIBDIR)/MUMPS_5.4.1
 LIB_MUMPS = -L$(MUMPS_DIR)/lib -ldmumps -lmumps_common -lpord
 INC_MUMPS = $(MUMPS_DIR)/include
 
 LIB_MUMPS += $(MKLLIB)
 INC_MUMPS_EXTRA += $(MKLINC)
 
endif
 
ifeq (1, $(USE_STRUMPACK))
 MKLLIB= -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 \
         -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64
 MKLINC = -I${MKLROOT}/include
 GKLIB_HOME = $(LIBDIR)/GKlib/install
 METIS_HOME = $(LIBDIR)/METIS/install
 PARMETIS_HOME = $(LIBDIR)/ParMETIS/install
 STRUMPACK_HOME=$(LIBDIR)/STRUMPACK/install
 
 ifeq (1, $(USE_INTSIZE64))
   METIS_HOME = $(LIB64DIR)/METIS_64/install
   PARMETIS_HOME = $(LIB64DIR)/ParMETIS_64/install
   STRUMPACK_HOME=$(LIB64DIR)/STRUMPACK_64/install
 endif
 
 STRUMPACKINC = -I$(STRUMPACK_HOME)/include
 STRUMPACKLIB = -L$(STRUMPACK_HOME)/lib64 -lstrumpack
 
 STRUMPACKINC += -I$(METIS_HOME)/include -I$(GKLIB_HOME)/include -I$(PARMETIS_HOME)/include
 STRUMPACKLIB += -L$(METIS_HOME)/lib -lmetis -L$(GKLIB_HOME)/lib -lGKlib -L$(PARMETIS_HOME)/lib -lparmetis
 
 STRUMPACKINC += $(MKLINC)
 STRUMPACKLIB += $(MKLLIB)
 
 DEFINES += -DNEWSPK
endif
 
# --- FFTW Library
INC_FFTW          = -I$(MKLROOT)/include/fftw
 
# --- HDF5 Library
HDF5INCLUDE = $(HDF5_HOME)/include
HDF5LIB     = -L$(HDF5_HOME)/lib -lhdf5_hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz
