// Parallel discrete cosine transform using MATAR + Kokkos.
//
// y:     real array of n elements (1-based DFMatrixKokkos, modified in-place).
// n:     number of elements (must be a power of 2).
// isign: +1 for forward transform, -1 for inverse.
//
// Parallelization:
//   Pre-transform pass: twiddle factors precomputed, DO_ALL over n/2-1
//   disjoint (j, n-j) pairs.  Sum contributions stored for later reduction.
//   Correction pass (isign==-1): parallel subtraction after serial sums.

#pragma once
#include <cmath>
#include <matar.h>
#include "realft.hpp"

using namespace mtr;

inline void cosft(DFMatrixKokkos<double>& y, int n, int isign)
{
    double theta = 3.14159265358979 / static_cast<double>(n);
    int m = n / 2;

    // Precompute twiddle factors for j=1..m-1
    int tw_len = (m - 1 > 0) ? m - 1 : 1;
    CArrayKokkos<double> tw_c(tw_len), tw_s(tw_len);
    CArrayKokkos<double> wr_y2(tw_len);

    if (m > 1) {
        FOR_ALL(j, 0, m - 1, {
            double angle = (j + 1) * theta;
            tw_c(j) = cos(angle);
            tw_s(j) = sin(angle);
        });
        MATAR_FENCE();
    }

    // Pre-transform: modify y pairs and store WR*Y2 contributions
    if (m > 1) {
        FOR_ALL(j, 0, m - 1, {
            double wi = tw_s(j);
            double wr = tw_c(j);
            int jp1 = j + 2;
            int nmj = n - j;
            double y1 = 0.5 * (y(jp1) + y(nmj));
            double y2 = y(jp1) - y(nmj);
            y(jp1) = y1 - wi * y2;
            y(nmj) = y1 + wi * y2;
            wr_y2(j) = wr * y2;
        });
        MATAR_FENCE();
    }

    // Compute sum = y(1) + sum(wr*y2) serially, store result
    CArrayKokkos<double> sum_store(1);
    int tw_count = m - 1;
    RUN({
        double s = y(1);
        for (int j = 0; j < tw_count; j++) {
            s += wr_y2(j);
        }
        sum_store(0) = s;
    });
    MATAR_FENCE();

    realft(y, m, 1);

    // Set y(2) = sum, then cumulative prefix sum on even indices
    RUN({
        double s = sum_store(0);
        y(2) = s;
        for (int j = 4; j <= n; j += 2) {
            s += y(j);
            y(j) = s;
        }
    });
    MATAR_FENCE();

    if (isign == -1) {
        // Compute EVEN and ODD sums, derive correction constants
        CArrayKokkos<double> params(2);
        RUN({
            double even_sum = y(1);
            double odd_sum  = y(2);
            for (int i = 3; i <= n - 1; i += 2) {
                even_sum += y(i);
                odd_sum  += y(i + 1);
            }
            double enf0 = 2.0 * (even_sum - odd_sum);
            double sumo = y(1) - enf0;
            double sume = (2.0 * odd_sum / n) - sumo;
            params(0) = sumo;
            params(1) = sume;
            y(1) = 0.5 * enf0;
            y(2) = y(2) - sume;
        });
        MATAR_FENCE();

        // Parallel correction of remaining elements
        int half = (n - 2) / 2;
        if (half > 0) {
            FOR_ALL(k, 0, half, {
                int i = 2 * k + 3;
                y(i)     -= params(0);
                y(i + 1) -= params(1);
            });
            MATAR_FENCE();
        }
    }
}
