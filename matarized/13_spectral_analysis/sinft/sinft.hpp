// Parallel discrete sine transform using MATAR + Kokkos.
//
// y: real array of n elements (1-based DFMatrixKokkos, modified in-place).
// n: number of elements (must be a power of 2).
//
// Parallelization:
//   Pre-transform pass: precomputed twiddle factors allow a DO_ALL
//   over n/2 disjoint (j, n-j) index pairs.
//   Prefix-sum pass: inherently serial (RUN on device).

#pragma once
#include <cmath>
#include <matar.h>
#include "realft.hpp"

using namespace mtr;

inline void sinft(DFMatrixKokkos<double>& y, int n)
{
    double theta = 3.14159265358979 / static_cast<double>(n);
    int m = n / 2;

    // Precompute sin(j*theta) for j=1..m
    CArrayKokkos<double> tw_s(m);
    FOR_ALL(j, 0, m, {
        tw_s(j) = sin((j + 1) * theta);
    });
    MATAR_FENCE();

    RUN({ y(1) = 0.0; });
    MATAR_FENCE();

    // Pre-transform: thread j touches y(j+2) and y(n-j), disjoint across threads
    FOR_ALL(j, 0, m, {
        double wi = tw_s(j);
        int jp1 = j + 2;
        int nmj = n - j;
        double y1 = wi * (y(jp1) + y(nmj));
        double y2 = 0.5 * (y(jp1) - y(nmj));
        y(jp1) = y1 + y2;
        y(nmj) = y1 - y2;
    });
    MATAR_FENCE();

    realft(y, m, 1);

    // Sequential prefix-sum rearrangement
    RUN({
        y(1) = 0.5 * y(1);
        y(2) = 0.0;
        double sum = 0.0;
        for (int j = 1; j <= n - 1; j += 2) {
            sum += y(j);
            y(j) = y(j + 1);
            y(j + 1) = sum;
        }
    });
    MATAR_FENCE();
}
