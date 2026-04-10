// Parallel evaluation of the maximum-entropy power spectrum
// using MATAR + Kokkos.
//
// fdt: frequency times delta-t (0 to 0.5).
// cof: LP coefficients from memcof, m elements (1-based DFMatrixKokkos).
// m:   number of LP coefficients.
// pm:  mean-square discrepancy from memcof.
//
// Returns: power spectral density at frequency fdt.
//
// Parallelization:
//   Twiddle factors precomputed via FOR_ALL.
//   Real/imag sums via two FOR_REDUCE_SUM reductions.

#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

inline double evlmem(double fdt, DFMatrixKokkos<double>& cof, int m, double pm)
{
    double theta = 6.28318530717959 * fdt;

    // Precompute twiddle factors cos(k*theta), sin(k*theta) for k=1..m
    CArrayKokkos<double> tw_r(m), tw_i(m);
    FOR_ALL(k, 0, m, {
        double angle = (k + 1) * theta;
        tw_r(k) = cos(angle);
        tw_i(k) = sin(angle);
    });
    MATAR_FENCE();

    // Parallel reductions for real and imaginary parts of denominator
    double sumr_part = 0.0, sumi_part = 0.0;
    double loc_r = 0.0, loc_i = 0.0;

    FOR_REDUCE_SUM(k, 0, m, loc_r, {
        loc_r += cof(k + 1) * tw_r(k);
    }, sumr_part);

    FOR_REDUCE_SUM(k, 0, m, loc_i, {
        loc_i += cof(k + 1) * tw_i(k);
    }, sumi_part);

    double sumr = 1.0 - sumr_part;
    double sumi = -sumi_part;

    return pm / (sumr * sumr + sumi * sumi);
}
