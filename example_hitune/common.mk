#
# common Makefile
#
# Gets included after the local Makefile in an example sub-directory
#
BASEDIR=$(dir $(lastword $(MAKEFILE_LIST)))
HYPERF_HOME ?= /root/hitune

# CC=gcc
#CC=icc
CC=clang

NPROCS=4
NTHREADS=4
POLYBENCH_UTIL_DIR ?= $(HYPERF_HOME)/polybench/OpenMP/utilities
POLYBENCHINCDIR := $(POLYBENCH_UTIL_DIR)
POLYBENCHSRC := $(POLYBENCH_UTIL_DIR)/polybench.c
PLC ?= $(HYPERF_HOME)/pluto/polycc
ifeq ("$(wildcard $(PLC))","")
PLC := $(BASEDIR)polycc
endif

# Intel MKL library paths
MKLROOT=/opt/intel/mkl

OPENBLAS_CFLAGS=-I/usr/include/openblas
OPENBLAS_LDFLAGS=-L/usr/lib64/openblas -lopenblas

BLIS_CFLAGS=-I/usr/include/blis
BLIS_LDFLAGS=-L/usr/local/lib/ -lblis

ifeq ($(CC), icc)
	OPT_FLAGS     := -O3 -xHost -ansi-alias -ipo -fp-model precise 
	PAR_FLAGS     := -parallel
	OMP_FLAGS     := -qopenmp
	MKL_CFLAGS    := -DMKL_ILP64 -mkl=parallel
	MKL_LDFLAGS   := -liomp5 -lpthread -lm -ldl
else ifeq ($(CC), clang)
	CLANG_RESOURCE_DIR := $(shell $(CC) -print-resource-dir 2>/dev/null)
	OMP_INCLUDE_DIR := $(CLANG_RESOURCE_DIR)/include
	OMP_LIB_DIR := /usr/local/lib/x86_64-unknown-linux-gnu
	OPT_FLAGS     := -I/usr/lib/ -I$(POLYBENCHINCDIR) -I. -I./../ -DPOLYBENCH_TIME -O3
	PAR_FLAGS     := --parallel 
	OMP_FLAGS     := -fopenmp=libomp -I$(OMP_INCLUDE_DIR) -L$(OMP_LIB_DIR) -Wl,-rpath,$(OMP_LIB_DIR) -lomp
	MKL_CFLAGS    := -DMKL_ILP64 -m64 -I$(MKLROOT)/include
	MKL_LDFLAGS   := -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
else
	# for gcc
	OPT_FLAGS     := -I$(POLYBENCHINCDIR) -I. -I./../ -DPOLYBENCH_TIME -O3 -march=skylake-avx512 
	PAR_FLAGS     := -ftree-parallelize-loops=$(NTHREADS) --parallel
	OMP_FLAGS     := -fopenmp 
	MKL_CFLAGS    := -DMKL_ILP64 -m64 -I$(MKLROOT)/include
	MKL_LDFLAGS   := -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl
endif

#CFLAGS += -DTIME
CFLAGS +=
LDFLAGS += -lm
PLCFLAGS += 
TILEFLAGS += --tile --smartfuse --prevector
BUILD_BIN ?= 0

ifdef POLYBENCH
	CFLAGS += -DPOLYBENCH_USE_SCALAR_LB -DPOLYBENCH_TIME -I $(POLYBENCHINCDIR) $(POLYBENCHSRC)
	DISTOPT_FLAGS += --variables_not_global
endif

all: par

$(SRC).opt.c:  $(SRC).c
	$(PLC) $(SRC).c --notile --noparallel $(PLCFLAGS)  -o $@

$(SRC).tiled.c:  $(SRC).c
	$(PLC) $(SRC).c --noparallel $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).par.c:  $(SRC).c
	$(PLC) $(SRC).c --parallel $(TILEFLAGS) $(PLCFLAGS)  -o $@

$(SRC).mlbpar.c:  $(SRC).c
	$(PLC) $(SRC).c --parallel --full-diamond-tile $(TILEFLAGS) $(PLCFLAGS)  -o $@

# Version that doesn't use diamond tiling
$(SRC).pipepar.c:  $(SRC).c
	$(PLC) $(SRC).c --parallel --nodiamond-tile $(TILEFLAGS) $(PLCFLAGS) -o $@

orig: $(SRC).c
ifeq ($(BUILD_BIN),1)
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).c $(POLYBENCHSRC) -o $@ $(LDFLAGS)
else
	@echo "[common.mk] BUILD_BIN=0: skip $@ binary compile (use BUILD_BIN=1)"
endif

orig_par: $(SRC).c
ifeq ($(BUILD_BIN),1)
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(PAR_FLAGS) $(SRC).c $(POLYBENCHSRC) -o $@ $(LDFLAGS)
else
	@echo "[common.mk] BUILD_BIN=0: skip $@ binary compile (use BUILD_BIN=1)"
endif

orig_omp: $(SRC).c
ifeq ($(BUILD_BIN),1)
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).c $(POLYBENCHSRC) -o $@ $(LDFLAGS)
else
	@echo "[common.mk] BUILD_BIN=0: skip $@ binary compile (use BUILD_BIN=1)"
endif

opt: $(SRC).opt.c
ifeq ($(BUILD_BIN),1)
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).opt.c $(POLYBENCHSRC) -o $@ $(LDFLAGS)
else
	@echo "[common.mk] BUILD_BIN=0: generated $(SRC).opt.c"
endif

tiled: $(SRC).tiled.c
ifeq ($(BUILD_BIN),1)
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(SRC).tiled.c $(POLYBENCHSRC) -o $@ $(LDFLAGS)
else
	@echo "[common.mk] BUILD_BIN=0: generated $(SRC).tiled.c"
endif

par: $(SRC).par.c
ifeq ($(BUILD_BIN),1)
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).par.c $(POLYBENCHSRC) -o $@  $(LDFLAGS)
else
	@echo "[common.mk] BUILD_BIN=0: generated $(SRC).par.c"
endif

mlbpar: $(SRC).mlbpar.c
ifeq ($(BUILD_BIN),1)
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).mlbpar.c $(POLYBENCHSRC) -o $@  $(LDFLAGS)
else
	@echo "[common.mk] BUILD_BIN=0: generated $(SRC).mlbpar.c"
endif

# Version that doesn't use diamond tiling
pipepar: $(SRC).pipepar.c
ifeq ($(BUILD_BIN),1)
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).pipepar.c -o $@  $(LDFLAGS)
else
	@echo "[common.mk] BUILD_BIN=0: generated $(SRC).pipepar.c"
endif

perf: orig tiled par orig_par
	rm -f .test
	./orig
	OMP_NUM_THREADS=$(NTHREADS) ./orig_par
	./tiled
	OMP_NUM_THREADS=$(NTHREADS) ./par 

# Compare performance with and without diamond tiling.
pipeperf: par pipepar
	rm -f .test
	OMP_NUM_THREADS=$(NTHREADS) ./par
	OMP_NUM_THREADS=$(NTHREADS) ./pipepar 

test: orig tiled par
	touch .test
	./orig 2> out_orig
	./tiled 2> out_tiled
	diff -q out_orig out_tiled
	OMP_NUM_THREADS=$(NTHREADS) ./par 2> out_par4
	rm -f .test
	diff -q out_orig out_par4
	@echo Success!

lbtest: par pipepar
	touch .test
	OMP_NUM_THREADS=$(NTHREADS) ./par 2> out_par4
	OMP_NUM_THREADS=$(NTHREADS) ./pipepar 2> out_pipepar4
	rm -f .test
	diff -q out_par4 out_pipepar4
	diff -q out_par4 out_fulldiamondtile4
	@echo Success!

opt-test: orig opt
	touch .test
	./orig > out_orig
	./opt > out_opt
	rm -f .test
	diff -q out_orig out_opt
	@echo Success!
	rm -f .test

clean:
	rm -f out_* *.pipepar.c *.tiled.c *.opt.c *.par.c orig opt tiled par sched orig_par \
		hopt hopt *.par2d.c *.out.* \
		*.kernel.* a.out $(EXTRA_CLEAN) tags tmp* gmon.out *~ .unroll \
	   	.vectorize par2d parsetab.py *.body.c *.pluto.c *.par.cloog *.tiled.cloog *.pluto.cloog

exec-clean:
	rm -f out_* opt orig tiled sched sched hopt hopt par pipepar orig_par *.out.* *.kernel.* a.out \
		$(EXTRA_CLEAN) tags tmp* gmon.out *~ par2d
