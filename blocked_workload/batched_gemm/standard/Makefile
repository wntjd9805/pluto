SRC= batched_gemm_standard_pluto

PLCFLAGS =  #--codegen-context=100
TILEFLAGS = #--second-level-tile

MKFILE_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
include $(MKFILE_DIR)../../common.mk

$(SRC).par2d.c: 
	$(PLC) $(SRC).c --tile --parallel --multipipe $(TILEFLAGS) $(PLCFLAGS)

par2d: $(SRC).par2d.c decls.h util.h
	$(CC) $(CFLAGS) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d

par2d_test: $(SRC).par2d.c decls.h util.h
	$(CC) $(OPT_FLAGS) -openmp -lm $(SRC).par2d.c -o par2d_test -DTEST

ptest: tiled_test par2d_test par_test
	./tiled_test 2> out_tiled
	$(OMP_RUN_ENV) ./par_test 2> out_par
	diff -q out_tiled out_par
	$(OMP_RUN_ENV) ./par2d_test 2> out_par2d
	diff -q out_tiled out_par2d
