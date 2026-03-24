#define N 32
#define BLOCK_SIZE_SMOOTH_DEFAULT 2
#define BLOCK_SIZE_RESIDUAL_DEFAULT 16
#define BLOCK_SIZE_RESTRICT_DEFAULT 16
#define BLOCK_SIZE_COARSE_DEFAULT 8
#define BLOCK_SIZE_PROLONG_DEFAULT 32

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif
#include "../../polybench.h"
#ifdef __cplusplus
}
#endif

#ifndef N
#error "N must be defined before including multigrid_shared.h"
#endif

#ifndef BLOCK_SIZE_SMOOTH_DEFAULT
#error "BLOCK_SIZE_SMOOTH_DEFAULT must be defined before including multigrid_shared.h"
#endif

#ifndef BLOCK_SIZE_RESIDUAL_DEFAULT
#error "BLOCK_SIZE_RESIDUAL_DEFAULT must be defined before including multigrid_shared.h"
#endif

#ifndef BLOCK_SIZE_RESTRICT_DEFAULT
#error "BLOCK_SIZE_RESTRICT_DEFAULT must be defined before including multigrid_shared.h"
#endif

#ifndef BLOCK_SIZE_COARSE_DEFAULT
#error "BLOCK_SIZE_COARSE_DEFAULT must be defined before including multigrid_shared.h"
#endif

#ifndef BLOCK_SIZE_PROLONG_DEFAULT
#error "BLOCK_SIZE_PROLONG_DEFAULT must be defined before including multigrid_shared.h"
#endif

#define NC (N / 2)
#define PRE_SMOOTH 2
#define POST_SMOOTH 2
#define COARSE_ITERS 10

#ifndef BLOCK_SIZE_SMOOTH
#define BLOCK_SIZE_SMOOTH BLOCK_SIZE_SMOOTH_DEFAULT
#endif
#ifndef BLOCK_SIZE_RESIDUAL
#define BLOCK_SIZE_RESIDUAL BLOCK_SIZE_RESIDUAL_DEFAULT
#endif
#ifndef BLOCK_SIZE_RESTRICT
#define BLOCK_SIZE_RESTRICT BLOCK_SIZE_RESTRICT_DEFAULT
#endif
#ifndef BLOCK_SIZE_COARSE
#define BLOCK_SIZE_COARSE BLOCK_SIZE_COARSE_DEFAULT
#endif
#ifndef BLOCK_SIZE_PROLONG
#define BLOCK_SIZE_PROLONG BLOCK_SIZE_PROLONG_DEFAULT
#endif

static void init_random_2d(double *A, int m, int n)
{
  for (int i = 0; i < m * n; ++i) {
    A[i] = (double)rand() / RAND_MAX - 0.5;
  }
}

static void init_zero_2d(double *A, int m, int n)
{
  for (int i = 0; i < m * n; ++i) {
    A[i] = 0.0;
  }
}

static void zero_dirichlet_boundary_fine(double A[N][N])
{
  for (int i = 0; i < N; ++i) {
    A[i][0] = 0.0;
    A[i][N - 1] = 0.0;
  }
  for (int j = 0; j < N; ++j) {
    A[0][j] = 0.0;
    A[N - 1][j] = 0.0;
  }
}

static void *xmalloc_align64(size_t bytes)
{
  void *ptr = NULL;
  int ret = posix_memalign(&ptr, 64, bytes);
  if (ret || !ptr) {
    fprintf(stderr, "alloc fail\n");
    exit(1);
  }
  return ptr;
}

static void smooth_gs_inplace(double u[N][N], double f[N][N], int sweeps)
{
#pragma scop
  for (int s = 0; s < sweeps; ++s) {
    for (int bi = 0; bi < (N + BLOCK_SIZE_SMOOTH - 1) / BLOCK_SIZE_SMOOTH; ++bi) {
      for (int bj = 0; bj < (N + BLOCK_SIZE_SMOOTH - 1) / BLOCK_SIZE_SMOOTH; ++bj) {
        for (int ii = 0; ii < BLOCK_SIZE_SMOOTH; ++ii) {
          for (int jj = 0; jj < BLOCK_SIZE_SMOOTH; ++jj) {
            if (bi * BLOCK_SIZE_SMOOTH + ii < N &&
                bj * BLOCK_SIZE_SMOOTH + jj < N) {
              if (bi * BLOCK_SIZE_SMOOTH + ii > 0 &&
                  bj * BLOCK_SIZE_SMOOTH + jj > 0 &&
                  bi * BLOCK_SIZE_SMOOTH + ii < N - 1 &&
                  bj * BLOCK_SIZE_SMOOTH + jj < N - 1) {
                u[bi * BLOCK_SIZE_SMOOTH + ii][bj * BLOCK_SIZE_SMOOTH + jj] =
                    (f[bi * BLOCK_SIZE_SMOOTH + ii][bj * BLOCK_SIZE_SMOOTH + jj] +
                     u[bi * BLOCK_SIZE_SMOOTH + ii - 1][bj * BLOCK_SIZE_SMOOTH + jj] +
                     u[bi * BLOCK_SIZE_SMOOTH + ii + 1][bj * BLOCK_SIZE_SMOOTH + jj] +
                     u[bi * BLOCK_SIZE_SMOOTH + ii][bj * BLOCK_SIZE_SMOOTH + jj - 1] +
                     u[bi * BLOCK_SIZE_SMOOTH + ii][bj * BLOCK_SIZE_SMOOTH + jj + 1]) *
                    0.25;
              }
            }
          }
        }
      }
    }
  }
#pragma endscop
}

static void compute_residual(double u[N][N], double f[N][N], double r[N][N])
{
  for (int i = 0; i < N; ++i) {
    r[i][0] = 0.0;
    r[i][N - 1] = 0.0;
  }
  for (int j = 0; j < N; ++j) {
    r[0][j] = 0.0;
    r[N - 1][j] = 0.0;
  }

    for (int bi = 0; bi < (N + BLOCK_SIZE_RESIDUAL - 1) / BLOCK_SIZE_RESIDUAL; ++bi) {
      for (int bj = 0; bj < (N + BLOCK_SIZE_RESIDUAL - 1) / BLOCK_SIZE_RESIDUAL; ++bj) {
#pragma scop
        for (int ii = 0; ii < BLOCK_SIZE_RESIDUAL; ++ii) {
          for (int jj = 0; jj < BLOCK_SIZE_RESIDUAL; ++jj) {
            if (bi * BLOCK_SIZE_RESIDUAL + ii > 0 &&
                bj * BLOCK_SIZE_RESIDUAL + jj > 0 &&
                bi * BLOCK_SIZE_RESIDUAL + ii < N - 1 &&
                bj * BLOCK_SIZE_RESIDUAL + jj < N - 1) {
              r[bi * BLOCK_SIZE_RESIDUAL + ii][bj * BLOCK_SIZE_RESIDUAL + jj] =
                  f[bi * BLOCK_SIZE_RESIDUAL + ii][bj * BLOCK_SIZE_RESIDUAL + jj] +
                  u[bi * BLOCK_SIZE_RESIDUAL + ii - 1][bj * BLOCK_SIZE_RESIDUAL + jj] +
                  u[bi * BLOCK_SIZE_RESIDUAL + ii + 1][bj * BLOCK_SIZE_RESIDUAL + jj] +
                  u[bi * BLOCK_SIZE_RESIDUAL + ii][bj * BLOCK_SIZE_RESIDUAL + jj - 1] +
                  u[bi * BLOCK_SIZE_RESIDUAL + ii][bj * BLOCK_SIZE_RESIDUAL + jj + 1] -
                  4.0 * u[bi * BLOCK_SIZE_RESIDUAL + ii][bj * BLOCK_SIZE_RESIDUAL + jj];
            }
          }
        }
#pragma endscop
      }
    }
}

static void restrict_full_weighting(double r_f[N][N], double f_c[NC][NC])
{
  for (int i = 0; i < NC; ++i) {
    f_c[i][0] = 0.0;
    f_c[i][NC - 1] = 0.0;
  }
  for (int j = 0; j < NC; ++j) {
    f_c[0][j] = 0.0;
    f_c[NC - 1][j] = 0.0;
  }

    for (int bI = 0; bI < (NC + BLOCK_SIZE_RESTRICT - 1) / BLOCK_SIZE_RESTRICT; ++bI) {
      for (int bJ = 0; bJ < (NC + BLOCK_SIZE_RESTRICT - 1) / BLOCK_SIZE_RESTRICT; ++bJ) {
#pragma scop
        for (int ii = 0; ii < BLOCK_SIZE_RESTRICT; ++ii) {
          for (int jj = 0; jj < BLOCK_SIZE_RESTRICT; ++jj) {
            if (bI * BLOCK_SIZE_RESTRICT + ii > 0 &&
                bJ * BLOCK_SIZE_RESTRICT + jj > 0 &&
                bI * BLOCK_SIZE_RESTRICT + ii < NC - 1 &&
                bJ * BLOCK_SIZE_RESTRICT + jj < NC - 1) {
              f_c[bI * BLOCK_SIZE_RESTRICT + ii][bJ * BLOCK_SIZE_RESTRICT + jj] =
                  (4.0 * r_f[2 * (bI * BLOCK_SIZE_RESTRICT + ii)][2 * (bJ * BLOCK_SIZE_RESTRICT + jj)] +
                   2.0 * (r_f[2 * (bI * BLOCK_SIZE_RESTRICT + ii) - 1][2 * (bJ * BLOCK_SIZE_RESTRICT + jj)] +
                          r_f[2 * (bI * BLOCK_SIZE_RESTRICT + ii) + 1][2 * (bJ * BLOCK_SIZE_RESTRICT + jj)] +
                          r_f[2 * (bI * BLOCK_SIZE_RESTRICT + ii)][2 * (bJ * BLOCK_SIZE_RESTRICT + jj) - 1] +
                          r_f[2 * (bI * BLOCK_SIZE_RESTRICT + ii)][2 * (bJ * BLOCK_SIZE_RESTRICT + jj) + 1]) +
                   (r_f[2 * (bI * BLOCK_SIZE_RESTRICT + ii) - 1][2 * (bJ * BLOCK_SIZE_RESTRICT + jj) - 1] +
                    r_f[2 * (bI * BLOCK_SIZE_RESTRICT + ii) - 1][2 * (bJ * BLOCK_SIZE_RESTRICT + jj) + 1] +
                    r_f[2 * (bI * BLOCK_SIZE_RESTRICT + ii) + 1][2 * (bJ * BLOCK_SIZE_RESTRICT + jj) - 1] +
                    r_f[2 * (bI * BLOCK_SIZE_RESTRICT + ii) + 1][2 * (bJ * BLOCK_SIZE_RESTRICT + jj) + 1])) *
                  (1.0 / 16.0);
            }
          }
        }
#pragma endscop
      }
    }
}

static void coarse_solve_gs_inplace(double e_c[NC][NC], double f_c[NC][NC], int iters)
{
#pragma scop
  for (int t = 0; t < iters; ++t) {
    for (int bi = 0; bi < (NC + BLOCK_SIZE_COARSE - 1) / BLOCK_SIZE_COARSE; ++bi) {
      for (int bj = 0; bj < (NC + BLOCK_SIZE_COARSE - 1) / BLOCK_SIZE_COARSE; ++bj) {
        for (int ii = 0; ii < BLOCK_SIZE_COARSE; ++ii) {
          for (int jj = 0; jj < BLOCK_SIZE_COARSE; ++jj) {
            if (bi * BLOCK_SIZE_COARSE + ii < NC &&
                bj * BLOCK_SIZE_COARSE + jj < NC) {
              if (bi * BLOCK_SIZE_COARSE + ii > 0 &&
                  bj * BLOCK_SIZE_COARSE + jj > 0 &&
                  bi * BLOCK_SIZE_COARSE + ii < NC - 1 &&
                  bj * BLOCK_SIZE_COARSE + jj < NC - 1) {
                e_c[bi * BLOCK_SIZE_COARSE + ii][bj * BLOCK_SIZE_COARSE + jj] =
                    (f_c[bi * BLOCK_SIZE_COARSE + ii][bj * BLOCK_SIZE_COARSE + jj] +
                     e_c[bi * BLOCK_SIZE_COARSE + ii - 1][bj * BLOCK_SIZE_COARSE + jj] +
                     e_c[bi * BLOCK_SIZE_COARSE + ii + 1][bj * BLOCK_SIZE_COARSE + jj] +
                     e_c[bi * BLOCK_SIZE_COARSE + ii][bj * BLOCK_SIZE_COARSE + jj - 1] +
                     e_c[bi * BLOCK_SIZE_COARSE + ii][bj * BLOCK_SIZE_COARSE + jj + 1]) *
                    0.25;
              }
            }
          }
        }
      }
    }
  }
#pragma endscop
}

static void prolong_and_correct(double u[N][N], double e_c[NC][NC])
{
    for (int bI = 0; bI < (NC + BLOCK_SIZE_PROLONG - 1) / BLOCK_SIZE_PROLONG; ++bI) {
      for (int bJ = 0; bJ < (NC + BLOCK_SIZE_PROLONG - 1) / BLOCK_SIZE_PROLONG; ++bJ) {
#pragma scop
        for (int ii = 0; ii < BLOCK_SIZE_PROLONG; ++ii) {
          for (int jj = 0; jj < BLOCK_SIZE_PROLONG; ++jj) {
            if (bI * BLOCK_SIZE_PROLONG + ii > 0 &&
                bJ * BLOCK_SIZE_PROLONG + jj > 0 &&
                bI * BLOCK_SIZE_PROLONG + ii < NC &&
                bJ * BLOCK_SIZE_PROLONG + jj < NC) {
              u[2 * (bI * BLOCK_SIZE_PROLONG + ii)][2 * (bJ * BLOCK_SIZE_PROLONG + jj)] +=
                  e_c[bI * BLOCK_SIZE_PROLONG + ii][bJ * BLOCK_SIZE_PROLONG + jj];
            }
          }
        }

        for (int ii = 0; ii < BLOCK_SIZE_PROLONG; ++ii) {
          for (int jj = 0; jj < BLOCK_SIZE_PROLONG; ++jj) {
            if (bI * BLOCK_SIZE_PROLONG + ii + 1 < NC &&
                bJ * BLOCK_SIZE_PROLONG + jj > 0 &&
                bJ * BLOCK_SIZE_PROLONG + jj < NC) {
              u[2 * (bI * BLOCK_SIZE_PROLONG + ii) + 1][2 * (bJ * BLOCK_SIZE_PROLONG + jj)] +=
                  0.5 * (e_c[bI * BLOCK_SIZE_PROLONG + ii][bJ * BLOCK_SIZE_PROLONG + jj] +
                         e_c[bI * BLOCK_SIZE_PROLONG + ii + 1][bJ * BLOCK_SIZE_PROLONG + jj]);
            }
          }
        }

        for (int ii = 0; ii < BLOCK_SIZE_PROLONG; ++ii) {
          for (int jj = 0; jj < BLOCK_SIZE_PROLONG; ++jj) {
            if (bJ * BLOCK_SIZE_PROLONG + jj + 1 < NC &&
                bI * BLOCK_SIZE_PROLONG + ii > 0 &&
                bI * BLOCK_SIZE_PROLONG + ii < NC) {
              u[2 * (bI * BLOCK_SIZE_PROLONG + ii)][2 * (bJ * BLOCK_SIZE_PROLONG + jj) + 1] +=
                  0.5 * (e_c[bI * BLOCK_SIZE_PROLONG + ii][bJ * BLOCK_SIZE_PROLONG + jj] +
                         e_c[bI * BLOCK_SIZE_PROLONG + ii][bJ * BLOCK_SIZE_PROLONG + jj + 1]);
            }
          }
        }

        for (int ii = 0; ii < BLOCK_SIZE_PROLONG; ++ii) {
          for (int jj = 0; jj < BLOCK_SIZE_PROLONG; ++jj) {
            if (bI * BLOCK_SIZE_PROLONG + ii + 1 < NC &&
                bJ * BLOCK_SIZE_PROLONG + jj + 1 < NC) {
              u[2 * (bI * BLOCK_SIZE_PROLONG + ii) + 1][2 * (bJ * BLOCK_SIZE_PROLONG + jj) + 1] +=
                  0.25 * (e_c[bI * BLOCK_SIZE_PROLONG + ii][bJ * BLOCK_SIZE_PROLONG + jj] +
                          e_c[bI * BLOCK_SIZE_PROLONG + ii + 1][bJ * BLOCK_SIZE_PROLONG + jj] +
                          e_c[bI * BLOCK_SIZE_PROLONG + ii][bJ * BLOCK_SIZE_PROLONG + jj + 1] +
                          e_c[bI * BLOCK_SIZE_PROLONG + ii + 1][bJ * BLOCK_SIZE_PROLONG + jj + 1]);
            }
          }
        }
#pragma endscop
      }
    }
}

static void vcycle_1level(
    double u[N][N], double f[N][N],
    double r[N][N],
    double e_c[NC][NC], double f_c[NC][NC])
{
  smooth_gs_inplace(u, f, PRE_SMOOTH);
  compute_residual(u, f, r);
  restrict_full_weighting(r, f_c);

#pragma scop
  for (int i = 0; i < NC; ++i) {
    for (int j = 0; j < NC; ++j) {
      e_c[i][j] = 0.0;
    }
  }
#pragma endscop

  coarse_solve_gs_inplace(e_c, f_c, COARSE_ITERS);
  prolong_and_correct(u, e_c);
  smooth_gs_inplace(u, f, POST_SMOOTH);
}

int main(void)
{
  srand(1234);

  double (*u)[N] = (double (*)[N])xmalloc_align64(sizeof(double) * N * N);
  double (*f)[N] = (double (*)[N])xmalloc_align64(sizeof(double) * N * N);
  double (*r)[N] = (double (*)[N])xmalloc_align64(sizeof(double) * N * N);
  double (*e_c)[NC] = (double (*)[NC])xmalloc_align64(sizeof(double) * NC * NC);
  double (*f_c)[NC] = (double (*)[NC])xmalloc_align64(sizeof(double) * NC * NC);

  init_zero_2d((double *)u, N, N);
  init_random_2d((double *)f, N, N);
  zero_dirichlet_boundary_fine(f);
  init_zero_2d((double *)r, N, N);
  init_zero_2d((double *)e_c, NC, NC);
  init_zero_2d((double *)f_c, NC, NC);

  vcycle_1level(u, f, r, e_c, f_c);

  polybench_timer_start();
  for (int it = 0; it < 15; ++it) {
    vcycle_1level(u, f, r, e_c, f_c);
  }
  polybench_timer_stop();
  polybench_timer_print();

  volatile double sink = u[1][1] + r[1][1] + f_c[1][1] + e_c[1][1];
  (void)sink;

  free(u);
  free(f);
  free(r);
  free(e_c);
  free(f_c);
  return 0;
}
