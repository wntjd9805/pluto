/**
 * jacobi-2d-imper.h: This file is part of the PolyBench/C 3.2 test suite.
 *
 *
 * Contact: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://polybench.sourceforge.net
 */
 #ifndef JACOBI_2D_IMPER_H
 # define JACOBI_2D_IMPER_H
 
 /* Default to STANDARD_DATASET. */
 # if !defined(MINI_DATASET) && !defined(SMALL_DATASET) && !defined(LARGE_DATASET) && !defined(EXTRALARGE_DATASET)
 #  define STANDARD_DATASET
 # endif
 
 /* Do not define anything if the user manually defines the size. */
 # if !defined(TSTEPS) && !defined(N)
 /* Define sample dataset sizes. */
 #  ifdef MINI_DATASET
 #   define TSTEPS 20
 #   define N 30
 #  endif
 
 #  ifdef SMALL_DATASET
 #   define TSTEPS 40
 #   define N 90
 #  endif
 
 #  ifdef STANDARD_DATASET
 #   define TSTEPS 100
 #   define N 250
 #  endif
 
 #  ifdef LARGE_DATASET
 #   define TSTEPS 500
 #   define N 1300
 #  endif
 
 #  ifdef EXTRALARGE_DATASET
 #   define TSTEPS 1000
 #   define N 2800
 #  endif
 
 # endif /* !N */
 
 # define _PB_TSTEPS POLYBENCH_LOOP_BOUND(TSTEPS,tsteps)
 # define _PB_N POLYBENCH_LOOP_BOUND(N,n)
 
 # ifndef DATA_TYPE
 #  define DATA_TYPE double
 #  define DATA_PRINTF_MODIFIER "%0.2lf "
 # endif
 
 
 #endif /* !JACOBI_2D_IMPER */