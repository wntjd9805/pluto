#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <string.h>

#define MAT_M   128
#define MAT_N2  128
#define KB      128

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 128
#endif

#ifdef __cplusplus
extern "C" {
#endif
#include "../../polybench.h"
#ifdef __cplusplus
}
#endif
/* ---------------------------------------------------------------------------
 * 64-byte aligned allocator
 * ------------------------------------------------------------------------- */
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


static void init_random(double *A, int m, int n)
{
    for (int i = 0; i < m * n; ++i)
        A[i] = (double)rand() / RAND_MAX;
}

/* ---------------------------------------------------------------------------
 * GEQRF trailing update
 * ------------------------------------------------------------------------- */
void geqrf_trailing_update_loop(
    double V[MAT_M][KB],
    double T[KB][KB],
    double C2[MAT_M][MAT_N2],
    double W1[KB][MAT_N2],
    double W2[KB][MAT_N2])
{
    // (1) W1 = V^T * C2
    for (int bi = 0; bi < (KB + BLOCK_SIZE - 1) / BLOCK_SIZE; ++bi) {
        for (int bk = 0; bk < (MAT_M + BLOCK_SIZE - 1) / BLOCK_SIZE; ++bk) {
            #pragma scop
            for (int ii = 0; ii < BLOCK_SIZE; ++ii) {
                for (int j = 0; j < MAT_N2; ++j) {
                    for (int kk = 0; kk < BLOCK_SIZE; ++kk) {

                        if ((bi*128 + ii) < KB &&
                            (bk*128 + kk) < MAT_M)
                        {
                            W1[bi*128 + ii][j] +=
                                V[bk*128 + kk][bi*128 + ii] *
                                C2[bk*128 + kk][j];
                        }
                    }
                }
            }
            #pragma endscop
        }
    }

     for (int bi = 0; bi < (KB + BLOCK_SIZE - 1) / BLOCK_SIZE; ++bi) 
        for (int bk = 0; bk < (KB + BLOCK_SIZE - 1) / BLOCK_SIZE; ++bk) {
            #pragma scop
            for (int ii = 0; ii < 128; ++ii) {
                for (int j = 0; j < MAT_N2; ++j) {
                    for (int kk = 0; kk < 128; ++kk) {
                        if ((bi*128 + ii) < KB &&
                            (bk*128 + kk) < KB)
                        {
                            W2[bi*128 + ii][j] +=
                                T[bi*128 + ii][bk*128 + kk] *
                                W1[bk*128 + kk][j];
                        }
                    }
                }
            }
            #pragma endscop
        }
     for (int bi = 0; bi < (MAT_M + BLOCK_SIZE - 1) / BLOCK_SIZE; ++bi) {
        for (int bk = 0; bk < (KB + BLOCK_SIZE - 1) / BLOCK_SIZE; ++bk) {
            #pragma scop
            for (int ii = 0; ii < 128; ++ii) {
                for (int j = 0; j < MAT_N2; ++j) {
                    for (int kk = 0; kk < 128; ++kk) {

                        if ((bi*128 + ii) < MAT_M &&
                            (bk*128 + kk) < KB)
                        {
                            C2[bi*128 + ii][j] -=
                                V[bi*128 + ii][bk*128 + kk] *
                                W2[bk*128 + kk][j];
                        }
                    }
                }
            }
            #pragma endscop
        }
    }
}

/* ---------------------------------------------------------------------------
 * Main
 * ------------------------------------------------------------------------- */
int main(void)
{
    srand(1234);

    // ============================================================
    // Allocate aligned matrices (64B aligned)
    // ============================================================
    double (*V)[KB]      = (double (*)[KB])xmalloc_align64(sizeof(double) * MAT_M * KB);
    double (*T)[KB]      = (double (*)[KB])xmalloc_align64(sizeof(double) * KB    * KB);
    double (*C2)[MAT_N2] = (double (*)[MAT_N2])xmalloc_align64(sizeof(double) * MAT_M * MAT_N2);

    double (*W1)[MAT_N2] = (double (*)[MAT_N2])xmalloc_align64(sizeof(double) * KB    * MAT_N2);
    double (*W2)[MAT_N2] = (double (*)[MAT_N2])xmalloc_align64(sizeof(double) * KB    * MAT_N2);

    init_random((double*)V, MAT_M, KB);
    init_random((double*)T, KB, KB);
    init_random((double*)C2, MAT_M, MAT_N2);

    // Benchmark
    geqrf_trailing_update_loop(V, T, C2, W1, W2);

    polybench_timer_start();
    for (int i = 0; i < 15; i++) {
        geqrf_trailing_update_loop(V, T, C2, W1, W2);
    }
    polybench_timer_stop();
    polybench_timer_print();

    // printf("Elapsed: %.6f sec\n", elapsed);
    // printf("C2[0,0] = %.6f, C2[last,last] = %.6f\n",
    //        C2[0][0], C2[MAT_M-1][MAT_N2-1]);

    // Free
    free(V);
    free(T);
    free(C2);
    free(W1);
    free(W2);

    return 0;
}