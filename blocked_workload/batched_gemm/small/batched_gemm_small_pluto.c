/******************************************************************************
 * File: batched_gemm_bench.c
 *
 * Purpose:
 *   Benchmark the **batched small GEMM** kernel—computing many tiny
 *   matrix multiplications in a batch-friendly loop:
 *
 *       For b = 0..batch-1:
 *         C_b = beta * C_b + alpha * A_b * B_b
 *
 *   This pattern appears in block factorizations on small tiles, BLAS batched
 *   APIs, and GPU/CPU libraries targeting high throughput for small matrices.
 *
 * Background:
 *   - When matrices are very small (e.g., 8–64), launching a full GEMM per op
 *     is inefficient. Grouping multiple small GEMMs into a **batched** loop
 *     improves locality and amortizes overhead (front-end, cache, TLB, etc.).
 *   - Many libraries (cuBLAS, MKL, OpenBLAS, BLIS, MAGMA) expose "batched"
 *     interfaces for this reason.
 *
 * What this benchmark does:
 *   - Allocates a batch of pointers to A_b (m×k), B_b (k×n), C_b (m×n)
 *   - Initializes them with random numbers (contiguous buffers with pointer
 *     arrays for low setup overhead)
 *   - Runs the batched GEMM kernel with a configurable **batch tile** size
 *   - Reports elapsed time (seconds)
 *
 * Parameters (macros):
 *   - MAT_BATCH:  number of GEMMs in the batch
 *   - MAT_M,N,K:  GEMM sizes (C: m×n, A: m×k, B: k×n)
 *   - BATCH_TILE: number of batch items processed per inner tile
 *   - ALPHA, BETA: scalar multipliers
 *
 * Build:
 *   gcc -O3 batched_gemm_bench.c -o batched_gemm_bench
 *
 * Run:
 *   ./batched_gemm_bench
 *
 * Notes:
 *   - Row-major layout is assumed with lda = K, ldb = N, ldc = N (matching
 *     the index expressions inside the kernel).
 *   - For real workloads, consider OpenMP parallelism across batch tiles or
 *     vectorization-friendly data layouts (AoSoA/interleave).
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
#define MAT_BATCH   1024 // number of small GEMMs in the batch
#define MAT_M       128    // rows of C and A
#define MAT_N       128    // cols of C and B
#define MAT_K       128    // cols of A / rows of B
// #define BATCH_TILE  32    // how many GEMMs to process per tile
#define ALPHA       1.0
#define BETA        0.0

#ifndef BATCH_TILE
#define BATCH_TILE 1024
#endif
// ---------------------------------------------------------------------------
// Kernel: multiple small (m×n, k) GEMMs in batch
//   Cp = beta·Cp + alpha·Ap·Bp
//   Row-major: Ap[i*lda + kk], Bp[kk*ldb + j], Cp[i*ldc + j]
// ---------------------------------------------------------------------------
void batched_small_gemm_loop(
     double A[MAT_BATCH][MAT_M][MAT_K],
     double B[MAT_BATCH][MAT_K][MAT_N],
    double       C[MAT_BATCH][MAT_M][MAT_N])

{   // Cp += alpha * Ap * Bp

if (BETA != 1.0) {

    for (int bi = 0; bi < (MAT_BATCH + BATCH_TILE - 1) / BATCH_TILE; ++bi) {
        #pragma scop
        for (int batch = 0; batch < 1024; ++batch) {
            for (int i = 0; i < MAT_M; ++i) {
                for (int j = 0; j < MAT_N; ++j) {
                    C[bi * 1024 + batch][i][j] *= BETA;
                        for (int kk = 0; kk < MAT_K; ++kk) { 
                            if (bi * 1024 + batch < MAT_BATCH) {
                            C[bi * 1024 + batch][i][j] += ALPHA * A[bi * 1024 + batch][i][kk] * B[bi * 1024 + batch][kk][j];
                        }
                    }
                }
            }
        }
        #pragma endscop
    }
} else {
    for (int bi = 0; bi < (MAT_BATCH + BATCH_TILE - 1) / BATCH_TILE; ++bi) {
        #pragma scop
        for (int batch = 0; batch < 1024; ++batch) {
            for (int i = 0; i < MAT_M; ++i) {
                for (int j = 0; j < MAT_N; ++j) {
                    for (int kk = 0; kk < MAT_K; ++kk) { 
                        if (bi * 1024 + batch < MAT_BATCH) {
                        C[bi * 1024 + batch][i][j] += ALPHA * A[bi * 1024 + batch][i][kk] * B[bi * 1024 + batch][kk][j];
                    }
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
static void init_random(double *M, size_t len)
{
    for (size_t i = 0; i < len; ++i)
        M[i] = (double)rand() / RAND_MAX;
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
    // srand((unsigned)time(NULL));
    srand(1234); 
    // printf("=== Batched Small GEMM Benchmark ===\n");
    // printf("Batch=%d, m×n×k=%d×%d×%d, tile=%d\n",
    //        MAT_BATCH, MAT_M, MAT_N, MAT_K, BATCH_TILE);

    const int m = MAT_M, n = MAT_N, k = MAT_K;
    const int lda = MAT_K;     // row-major: A(m×k)
    const int ldb = MAT_N;     // row-major: B(k×n)
    const int ldc = MAT_N;     // row-major: C(m×n)
    const size_t elements_A = (size_t)MAT_BATCH * MAT_M * MAT_K;
    const size_t elements_B = (size_t)MAT_BATCH * MAT_K * MAT_N;
    const size_t elements_C = (size_t)MAT_BATCH * MAT_M * MAT_N;
    double (*A)[MAT_M][MAT_K] =
        (double (*)[MAT_M][MAT_K])xmalloc_align64(sizeof(double) * elements_A);

    double (*B)[MAT_K][MAT_N] =
        (double (*)[MAT_K][MAT_N])xmalloc_align64(sizeof(double) * elements_B);

    double (*C)[MAT_M][MAT_N] =
        (double (*)[MAT_M][MAT_N])xmalloc_align64(sizeof(double) * elements_C);
        

    if (!A || !B || !C) {
        fprintf(stderr, "Memory allocation failed (buffers).\n");
        free(A); free(B); free(C);
        return -1;
    }

    // Pointer arrays to each per-batch matrix
    // Initialize data
    init_random((double *)A, elements_A);
    init_random((double *)B, elements_B);
    init_random((double *)C, elements_C);
    batched_small_gemm_loop(
        A,
        B,
        C
    );
    polybench_timer_start();
    for (int i = 0; i < 15; i++) {
        batched_small_gemm_loop(
        A,
        B,
        C
    );
    }
    polybench_timer_stop();
    polybench_timer_print();

    // printf("Elapsed time: %.6f seconds\n", sec);

    // Print a couple of sample outputs for a sanity check
    double c00 = C[0][0][0];
    double clast = C[MAT_BATCH-1][m-1][n-1];
    // printf("Sample C[0](0,0)=%.6f, C[last](m-1,n-1)=%.6f\n", c00, clast);

    free(A); free(B); free(C);
    return 0;
}
