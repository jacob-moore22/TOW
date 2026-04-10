// Parallel maximum-entropy method (Burg) coefficient estimation
// using MATAR + Kokkos.
//
// data:  input signal of n elements (1-based DFMatrixKokkos).
// n:     signal length.
// m:     number of LP coefficients to compute.
// pm:    output mean-square discrepancy (scalar, host).
// cof:   output LP coefficients, m elements (1-based DFMatrixKokkos).
// wk1, wk2: work arrays of n elements each.
// wkm:  work array of m elements.
//
// Parallelization:
//   Initial power -- DO_REDUCE_SUM over n elements.
//   PNEUM/DENOM  -- DO_REDUCE_SUM over n-k elements per stage.
//   WKM copy     -- DO_ALL.
//   WK1/WK2 update -- serial (cross-element dependency).

#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

inline void memcof(DFMatrixKokkos<double>& data, int n, int m,
                   double& pm,
                   DFMatrixKokkos<double>& cof,
                   DFMatrixKokkos<double>& wk1,
                   DFMatrixKokkos<double>& wk2,
                   DFMatrixKokkos<double>& wkm)
{
    // Initial power estimate
    double p_sum = 0.0;
    double loc_p = 0.0;
    DO_REDUCE_SUM(j, 1, n, loc_p, {
        loc_p += data(j) * data(j);
    }, p_sum);
    pm = p_sum / n;

    // Initialize work arrays
    RUN({
        wk1(1) = data(1);
        wk2(n - 1) = data(n);
    });
    MATAR_FENCE();

    if (n > 2) {
        DO_ALL(j, 2, n - 1, {
            wk1(j) = data(j);
            wk2(j - 1) = data(j);
        });
        MATAR_FENCE();
    }

    for (int k = 1; k <= m; k++) {
        // Parallel reductions for numerator and denominator
        double pneum = 0.0, denom = 0.0;
        double loc_num = 0.0, loc_den = 0.0;
        int nk = n - k;

        DO_REDUCE_SUM(j, 1, nk, loc_num, {
            loc_num += wk1(j) * wk2(j);
        }, pneum);

        DO_REDUCE_SUM(j, 1, nk, loc_den, {
            loc_den += wk1(j) * wk1(j) + wk2(j) * wk2(j);
        }, denom);

        double cof_k = 2.0 * pneum / denom;
        pm *= (1.0 - cof_k * cof_k);

        // Update cof(1..k-1) from previous wkm values
        if (k > 1) {
            int km1 = k - 1;
            RUN({
                for (int i = 1; i <= km1; i++)
                    cof(i) = wkm(i) - cof_k * wkm(km1 + 1 - i);
            });
            MATAR_FENCE();
        }

        RUN({ cof(k) = cof_k; });
        MATAR_FENCE();

        if (k == m) return;

        // Copy cof to wkm
        DO_ALL(i, 1, k, { wkm(i) = cof(i); });
        MATAR_FENCE();

        // Update wk1, wk2 (serial: cross-element read dependency)
        int nk1 = n - k - 1;
        RUN({
            double wmk = wkm(k);
            for (int j = 1; j <= nk1; j++) {
                wk1(j) = wk1(j) - wmk * wk2(j);
                wk2(j) = wk2(j + 1) - wmk * wk1(j + 1);
            }
        });
        MATAR_FENCE();
    }
}
