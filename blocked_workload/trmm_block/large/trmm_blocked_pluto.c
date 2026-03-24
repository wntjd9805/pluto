/* POLYBENCH/GPU-OPENMP
 *
 * This file is a part of the Polybench/GPU-OpenMP suite
 *
 * Contact:
 * William Killian <killian@udel.edu>
 *
 * Copyright 2013, The University of Delaware
 */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
/* Include polybench common header. */
#include "../../polybench.h"

/* Include benchmark-specific header. */
/* Default data type is double, default size is 4000. */
#include "../../block.h"

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef BLOCK_SIZE
#  define BLOCK_SIZE 8
#endif

/* Array initialization. */
static void init_array(int ni, DATA_TYPE *alpha,
                       DATA_TYPE POLYBENCH_2D(A, NI, NI, ni, ni),
                       DATA_TYPE POLYBENCH_2D(B, NI, NI, ni, ni)) {
  int i, j;

  *alpha = 32412;
  for (i = 0; i < ni; i++)
    for (j = 0; j < ni; j++) {
      A[i][j] = ((DATA_TYPE)i * j) / ni;
      B[i][j] = ((DATA_TYPE)i * j) / ni;
    }
}

/* DCE code. Must scan the entire live-out data.
  Can be used also to check the correctness of the output. */
static void print_array(int ni, DATA_TYPE POLYBENCH_2D(B, NI, NI, ni, ni)) {
  int i, j;

  for (i = 0; i < ni; i++)
    for (j = 0; j < ni; j++) {
      fprintf(stderr, DATA_PRINTF_MODIFIER, B[i][j]);
      if ((i * ni + j) % 20 == 0)
        fprintf(stderr, "\n");
    }
  fprintf(stderr, "\n");
}


static void kernel_trmm_blocked(int ni, DATA_TYPE alpha,
                                DATA_TYPE POLYBENCH_2D(A, NI, NI, ni, ni),
                                DATA_TYPE POLYBENCH_2D(B, NI, NI, ni, ni)) {
  const int n = ni;
  size_t alignment = 64;
  size_t size = (size_t)_PB_NI * _PB_NI * sizeof(DATA_TYPE);

  // 사이즈를 alignment의 배수로 맞춤 (안전을 위해)
  size_t aligned_size = (size + alignment - 1) & ~(alignment - 1);

  DATA_TYPE *B_read = (DATA_TYPE *)aligned_alloc(alignment, aligned_size);

    {  
    #pragma scop
    for (int i = 0; i < _PB_NI; ++i){
      for (int j = 0; j < _PB_NI; ++j)
        B_read[i * 2048 + j] = B[i][j];
    }
    #pragma endscop
    }

    for (int bi = 0; bi < (_PB_NI + BLOCK_SIZE - 1) / BLOCK_SIZE; ++bi) {
      for (int bk = bi + 1; bk < (_PB_NI + BLOCK_SIZE - 1) / BLOCK_SIZE; ++bk) {
        {
        #pragma scop
        for (int ii = 0; ii < BLOCK_SIZE; ++ii) {
          for (int j = 0; j < _PB_NI; ++j) {
            for (int kk = 0; kk < BLOCK_SIZE; ++kk) {
              if ((bi * 8 + ii) < _PB_NI &&
                  (bk * 8 + kk) < _PB_NI) {
                B[bi * 8 + ii][j] +=
                    A[bk * 8 + kk][bi * 8 + ii] *
                    B_read[(bk * 8 + kk) * 2048 + j];
              }
            }
          }
        }
        #pragma endscop
        }
      }

      {
      #pragma scop
      for (int ii = 0; ii < BLOCK_SIZE; ++ii) {
        for (int j = 0; j < _PB_NI; ++j) {
          for (int t = ii + 1; t < BLOCK_SIZE; ++t) {
            if ((bi * 8 + ii) < _PB_NI && (bi * 8 + t) < _PB_NI) {
              B[bi * 8 + ii][j] +=
                  A[bi * 8 + t][bi * 8 + ii] *
                  B_read[(bi * 8 + t) * 2048 +j]; // 읽기는 스냅샷, 쓰기는 B
            }
          }
        }
      }
      #pragma endscop
    }
    }

    {
    #pragma scop
    for (int i = 0; i < _PB_NI; ++i)
      for (int j = 0; j < _PB_NI; ++j)
        B[i][j] *= alpha;
    #pragma endscop
    }

  free(B_read);
}

int main(int argc, char **argv) {
  /* Retrieve problem size. */
  int ni = NI;

  /* Variable declaration/allocation. */
  DATA_TYPE alpha;
  POLYBENCH_2D_ARRAY_DECL(A,DATA_TYPE,NI,NI,ni,ni);
  POLYBENCH_2D_ARRAY_DECL(B,DATA_TYPE,NI,NI,ni,ni);

  /* Initialize array(s). */
  init_array(ni, &alpha, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  kernel_trmm_blocked(ni, alpha, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Run kernel. */
  polybench_timer_start();
  for (int i = 0; i < 15; i++) {
    kernel_trmm_blocked(ni, alpha, POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));
  }
  polybench_timer_stop();

  polybench_timer_print();
  // print_array(10, POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
    by the function call in argument. */
  polybench_prevent_dce(print_array(ni, POLYBENCH_ARRAY(B)));

  /* Be clean. */
  // POLYBENCH_FREE_ARRAY(A);
  // POLYBENCH_FREE_ARRAY(B);

  return 0;
}