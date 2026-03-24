/******************************************************************************
 * File: conv_im2col_gemm_blocked.c
 *
 * Purpose:
 *   Benchmark a forward convolution implemented using the classic
 *   im2col + GEMM strategy, but with blocking (tiling) over the
 *   output spatial dimension H_out * W_out.
 *
 *   All tensors are stored as fixed-size 2D or 3D arrays, defined by
 *   compile-time macros (true static arrays, no dynamic allocation).
 *
 * Blocking strategy:
 *   Let HW = H_out * W_out. We choose a tile size TILE_HW and process
 *   the output spatial positions in chunks:
 *
 *     for (hw0 = 0; hw0 < HW; hw0 += TILE_HW) {
 *       hw_len = min(TILE_HW, HW - hw0);
 *       im2col_nchw_tile(hw0, hw_len);   // fill col_tile[KIC][hw_len]
 *       gemm_nn_tile(hw0, hw_len);       // out_[oc][hw] += Wpk * col_tile
 *     }
 *
 *   This way, we only materialize a small im2col buffer for each tile
 *   instead of the entire Col[KIC][HW] at once, and GEMM is also done
 *   tile-by-tile for better cache behavior.
 *
 * Build:
 *   gcc -O3 conv_im2col_gemm_blocked.c -o conv_bench
 *
 * Run:
 *   ./conv_bench
 *
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
/* ---------------------------------------------------------------------------
 * Problem configuration
 * ------------------------------------------------------------------------- */
#define IN_C   64
#define IN_H   112
#define IN_W   112

#define OUT_C  128
#define K_H     3
#define K_W     3
#define STRIDE  1
#define PAD     1

#define H_OUT  112
#define W_OUT  112
#define HW_OUT 12544

/* Inner dimension for GEMM: K * C_in */
// 64*3*3=576 
#define KIC  576

/* Tile size for HW dimension (H_out * W_out) */
#ifndef TILE_HW
#define TILE_HW 4096
#endif
/* ---------------------------------------------------------------------------
 * Arrays (true 2D and 3D arrays, NOT dynamic)
 * ------------------------------------------------------------------------- */

/* Input tensor: x[C_in][H][W] */
static double x[IN_C][IN_H][IN_W];

/* Weights: w[C_out][C_in][K_h][K_w] */
static double w[OUT_C][IN_C][K_H][K_W];

/* Packed weights: Wpk[C_out][KIC] */
static double Wpk[OUT_C][KIC];

/*
 * Tile buffer for im2col:
 *   For each tile, we store col_tile[KIC][hw_len], where hw_len <= TILE_HW.
 *   We allocate the full [KIC][TILE_HW] statically and use only hw_len columns.
 */
static double col_tile[KIC][8192];

// #((112 + 2*1 - 3)/1 + 1) * ((112 + 2*1 - 3)/1 + 1) = 12544

/* Output matrix: out_[C_out][H_out * W_out] */
static double out_[OUT_C][12544]; /* named out_ to avoid 'out' keyword */


/* ---------------------------------------------------------------------------
 * Utility: Random initialization
 * ------------------------------------------------------------------------- */
static void random_fill(double *p, size_t n)
{
    for (size_t i = 0; i < n; ++i)
        p[i] = (double)rand() / (double)RAND_MAX - 0.5;
}

/* ---------------------------------------------------------------------------
 * im2col for a tile of HW
 *
 *  Input:
 *    x[IN_C][IN_H][IN_W]
 *
 *  Output:
 *    col_tile[KIC][hw_len]   (row-major, but stored as 2D static array)
 *
 *  We process output positions with global index:
 *    hw = oh * W_OUT + ow,  hw in [hw_start, hw_start + hw_len)
 * ------------------------------------------------------------------------- */
static void im2col(void)

{// 1. 바깥 루프는 타일 개수만큼 반복 (Stride 1)
    for (int oh = 0; oh < 112; ++oh) {
        for (int ow0_idx = 0; ow0_idx < (112 + 4096 - 1) / 4096; ++ow0_idx) {
#pragma scop
            for (int idx = 0; idx < 4096; ++idx) {
                    for (int c = 0; c < IN_C; ++c) {
                        for (int kh = 0; kh < K_H; ++kh) {
                            for (int kw = 0; kw < K_W; ++kw) {
                                if (ow0_idx * 4096 + idx < 112) {
                                if (0 <= oh - 1 + kh && oh - 1 + kh < 112 &&
                                    0 <= ow0_idx * 4096 + idx - 1 + kw &&
                                    ow0_idx * 4096 + idx - 1 + kw < 112) {
                                    col_tile[c * 9 + kh * 3 + kw][idx] =
                                        x[c][oh - 1 + kh][ow0_idx * 4096 + idx - 1 + kw];
                                } else {
                                    col_tile[c * 9 + kh * 3 + kw][idx] = 0;
                                }
                            }
                        }
                    }
                }
            }
            for (int oc = 0; oc < OUT_C; ++oc) {
                for (int idx = 0; idx < 4096; ++idx) {
                    if (ow0_idx * 4096 + idx < 112) {
                        for (int k = 0; k < KIC; ++k) {
                            out_[oc][oh * 112 + ow0_idx * 4096 + idx] +=
                                Wpk[oc][k] * col_tile[k][idx];
                        }
                    }
                }
            }
#pragma endscop
        }
    }
}

/* ---------------------------------------------------------------------------
 * Pack weights into 2D matrix (Wpk[OUT_C][KIC])
 * ------------------------------------------------------------------------- */
static void pack_weights(void)
{
    const int K = K_H * K_W;

    for (int oc = 0; oc < OUT_C; ++oc) {
        for (int ic = 0; ic < IN_C; ++ic) {
            for (int kh = 0; kh < K_H; ++kh) {
                for (int kw = 0; kw < K_W; ++kw) {

                    Wpk[oc][ic*K + kh*K_W + kw] = w[oc][ic][kh][kw];
                }
            }
        }
    }
}

/* ---------------------------------------------------------------------------
 * MAIN
 * ------------------------------------------------------------------------- */
int main(void)
{
    // srand((unsigned)time(NULL));
    srand(1234);

    // printf("=== Conv2D via im2col + GEMM – Blocked Static 2D Arrays Version ===\n");
    // printf("Input : %d x %d x %d\n", IN_C, IN_H, IN_W);
    // printf("Kernel: %d x %d, OUT_C=%d\n", K_H, K_W, OUT_C);
    // printf("Output: H_out=%d, W_out=%d (HW=%d)\n", H_OUT, W_OUT, HW_OUT);
    // printf("KIC   : %d\n", KIC);
    // printf("TILE_HW: %d\n\n", TILE_HW);

    /* Initialize arrays */
    random_fill(&x[0][0][0], (size_t)IN_C * IN_H * IN_W);
    random_fill(&w[0][0][0][0], (size_t)OUT_C * IN_C * K_H * K_W);

    /* Pack weights (not timed separately, but you can time if you want) */
    pack_weights();

    /* Benchmark im2col + GEMM with HW blocking */
    im2col();
    polybench_timer_start();
    for (int i = 0; i < 15; i++) {
        im2col();
    }
    polybench_timer_stop();
    polybench_timer_print();

    // printf("Sample out_[0][0]      = %.6f\n", out_[0][0]);
    // printf("Sample out_[last][last]= %.6f\n",
    //        out_[OUT_C-1][H_OUT*W_OUT - 1]);

    return 0;
}