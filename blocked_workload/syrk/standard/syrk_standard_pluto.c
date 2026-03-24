/******************************************************************************
 * File: syrk_bench.c
 *
 * Purpose:
 *   This program benchmarks the lower-triangular **SYRK**-like trailing update:
 *       C(lower) ← C(lower) - A * Aᵀ
 *   where A is (N2 × KB) and only the lower triangle of C is updated.
 *
 * Background:
 *   This kernel is the core of many blocked factorizations and updates:
 *     - Cholesky (POTRF) trailing update (A22 -= L21 * L21ᵀ)
 *     - Covariance / Gram matrix updates
 *     - Schur complement updates in various decompositions
 *
 *   It is equivalent to a symmetric rank-k update (BLAS-3 SYRK), and typically
 *   dominates runtime for large problems, thus benefiting from cache-blocking
 *   and GEMM-level optimizations in BLAS libraries (OpenBLAS, MKL, BLIS, MAGMA).
 *
 * What this benchmark does:
 *   - Randomly initializes A21 (N2 × KB) and C22 (N2 × N2)
 *   - Runs the lower-triangular SYRK update:
 *         for i = 0..N2-1
 *           for j = 0..i
 *             C22[i,j] -= Σ_p A21[i,p] * A21[j, p]
 *   - Measures and prints elapsed execution time
 *
 * Parameters (macros):
 *   - MAT_N2:  number of rows/cols in C22 (C22: MAT_N2 × MAT_N2)
 *   - KB:      rank parameter (A21: MAT_N2 × KB)
 *
 * Build:
 *   gcc -O3 syrk_bench.c -o syrk_bench
 *
 * Run:
 *   ./syrk_bench
 *
 * Author: Oh-kyoung Kwon
 * Date:   2025-10-30
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#ifdef __cplusplus
extern "C" {
#endif
#include "../../polybench.h"
#ifdef __cplusplus
}
#endif

#define MAT_N2  1024   // Size of C22 (C22: MAT_N2 × MAT_N2)
#define KB       1024   // Rank-k (A21: MAT_N2 × KB)

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 128
#endif

// ---------------------------------------------------------------------------
// Kernel function: SYRK lower update
//   C(lower) ← C(lower) - A * Aᵀ
//   A: N2 × kb, C: N2 × N2 (only lower triangle updated)
// ---------------------------------------------------------------------------
void syrk_lower_update_loop(
    double A21[MAT_N2][KB],
    double       C22[MAT_N2][MAT_N2])
{
        for (int bi = 0; bi < (MAT_N2 + BLOCK_SIZE - 1) / BLOCK_SIZE; ++bi) {
            for (int bj = 0; bj <= bi; ++bj) {
                #pragma scop
                for (int ii = 0; ii < BLOCK_SIZE; ++ii) {
                    for (int jj = 0; jj < BLOCK_SIZE; ++jj) {
                        for (int kk = 0; kk < KB; ++kk) {
                            if (bi * 128 + ii < MAT_N2 &&
                                bj * 128 + jj < MAT_N2 &&
                                bj * 128 + jj <= bi * 128 + ii)
                            {
                                C22[bi * 128 + ii][bj * 128 + jj] -=
                                    A21[bi * 128 + ii][kk] *
                                    A21[bj * 128 + jj][kk];
                            }
                        }
                    }
                }
                #pragma endscop
            }
        }
        
}

// ---------------------------------------------------------------------------
// Utilities
// ---------------------------------------------------------------------------
static void init_random(double *A, int m, int n)
{
    for (int i = 0; i < m * n; ++i)
        A[i] = (double)rand() / RAND_MAX;
}

static void *xmalloc_align64(size_t bytes)
{
    void *ptr = NULL;
    int ret = posix_memalign(&ptr, 64, bytes);
    if (ret || !ptr) {
        fprintf(stderr, "64-byte aligned allocation failed\n");
        exit(1);
    }
    return ptr;
}


// ---------------------------------------------------------------------------
// Main benchmark driver
// ---------------------------------------------------------------------------
int main(void)
{
    // srand(time(NULL));
    srand(1234); 
    // printf("=== SYRK Lower Update Benchmark ===\n");
    // printf("A21: %d×%d,  C22(lower): %d×%d\n", MAT_N2, KB, MAT_N2, MAT_N2);

    // Row-major layout with lda = number of columns for simplicity here.
    // We treat lda/ldc as leading dimensions consistent with the loops above.
    double (*A21)[KB] = (double (*)[KB])xmalloc_align64(sizeof(double) * MAT_N2 * KB);
    double (*C22)[MAT_N2] = (double (*)[MAT_N2])xmalloc_align64(sizeof(double) * MAT_N2 * MAT_N2);


    init_random((double*)A21, MAT_N2, KB);
    init_random((double*)C22, MAT_N2, MAT_N2);

    syrk_lower_update_loop(A21, C22);

    polybench_timer_start();
    for (int i = 0; i < 15; i++) {
        syrk_lower_update_loop(A21, C22);
    }
    polybench_timer_stop();
    polybench_timer_print();

    // printf("Sample C22[0,0] = %.6f, C22[last,last] = %.6f\n",
    //        C22[0][0], C22[MAT_N2-1][MAT_N2-1]);
    free(A21);
    free(C22);

    return 0;
}
