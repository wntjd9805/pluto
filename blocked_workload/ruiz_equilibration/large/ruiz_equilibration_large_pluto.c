/******************************************************************************
 * File: ruiz_equilibration_bench.c
 *
 * Fused norms (row_sq + col_sq in one pass) + compute r,c + apply scaling.
 *
 * Build:
 *   clang -O3 -fopenmp ruiz_equilibration_bench.c -o ruiz_bench -lm
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif
#include "../../polybench.h"
#ifdef __cplusplus
}
#endif

#define MAT_N   2048
#define EPS_VAL 1e-12

#ifndef BLOCK_SIZE_NORM
#define BLOCK_SIZE_NORM 2048
#endif

#ifndef BLOCK_SIZE_SCALE
#define BLOCK_SIZE_SCALE 2048
#endif

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------
static void init_random(double *A, int m, int n) {
  for (int i = 0; i < m * n; ++i) A[i] = (double)rand() / RAND_MAX - 0.5;
}

static void *xmalloc_align64(size_t bytes) {
  void *ptr = NULL;
  int ret = posix_memalign(&ptr, 64, bytes);
  if (ret || !ptr) {
    fprintf(stderr, "64-byte aligned allocation failed\n");
    exit(1);
  }
  return ptr;
}

// ---------------------------------------------------------------------------
// Kernel: norms + scaling (single equilibration step)
// ---------------------------------------------------------------------------
void ruiz_equilibration(
    double A[MAT_N][MAT_N],
    double row_sq[MAT_N],
    double col_sq[MAT_N],
    double r[MAT_N],
    double c[MAT_N])
{
  {
    #pragma scop
    // (0) init row_sq / col_sq
    for (int i = 0; i < MAT_N; ++i) {
      row_sq[i] = 0.0;
      col_sq[i] = 0.0;
    }

    // (1) norms: one pass over A, blocked
    for (int bi = 0; bi < 1; ++bi) {
      for (int bj = 0; bj < 1; ++bj) {
        for (int ii = 0; ii < 2048; ++ii) {
          for (int jj = 0; jj < 2048; ++jj) {
            // 가장 inner에서 범위 체크
            if (bi * 2048 + ii < MAT_N &&
                bj * 2048 + jj < MAT_N) {
              row_sq[bi * 2048 + ii] += A[bi * 2048 + ii][bj * 2048 + jj] * A[bi * 2048 + ii][bj * 2048 + jj];
              col_sq[bj * 2048 + jj] += A[bi * 2048 + ii][bj * 2048 + jj] * A[bi * 2048 + ii][bj * 2048 + jj];
            }
          }
        }
      }
    }

    // (2) compute r,c
    for (int i = 0; i < MAT_N; ++i) {
      r[i] = 1.0 / sqrt(row_sq[i] + EPS_VAL);
      c[i] = 1.0 / sqrt(col_sq[i] + EPS_VAL);
    }

    // (3) apply scaling in-place (blocked)
    for (int bi = 0; bi < 1; ++bi) {
      for (int bj = 0; bj < 1; ++bj) {
        for (int ii = 0; ii < 2048; ++ii) {
          for (int jj = 0; jj < 2048; ++jj) {
            // 가장 inner에서 범위 체크
            if (bi * 2048 + ii < MAT_N &&
                bj * 2048 + jj < MAT_N) {
              A[bi * 2048 + ii][bj * 2048 + jj] =
                  r[bi * 2048 + ii] *
                  A[bi * 2048 + ii][bj * 2048 + jj] *
                  c[bj * 2048 + jj];
            }
          }
        }
      }
    }
    #pragma endscop
  }
}

// ---------------------------------------------------------------------------
// Main benchmark driver
// ---------------------------------------------------------------------------
int main(void)
{
  srand(1234);

  double (*A)[MAT_N] = (double (*)[MAT_N])xmalloc_align64(sizeof(double) * MAT_N * MAT_N);
  double *row_sq     = (double *)xmalloc_align64(sizeof(double) * MAT_N);
  double *col_sq     = (double *)xmalloc_align64(sizeof(double) * MAT_N);
  double *r          = (double *)xmalloc_align64(sizeof(double) * MAT_N);
  double *c          = (double *)xmalloc_align64(sizeof(double) * MAT_N);

  init_random((double*)A, MAT_N, MAT_N);

  // warm-up
  ruiz_equilibration(A, row_sq, col_sq, r, c);

  polybench_timer_start();
  for (int it = 0; it < 15; ++it) {
    ruiz_equilibration(A, row_sq, col_sq, r, c);
  }
  polybench_timer_stop();
  polybench_timer_print();

  volatile double sink = A[0][0] + row_sq[0] + col_sq[0] + r[0] + c[0];
  (void)sink;

  free(A);
  free(row_sq);
  free(col_sq);
  free(r);
  free(c);
  return 0;
}
