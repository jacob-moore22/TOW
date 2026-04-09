// Parallel FFT of a single real-valued sequence using MATAR + Kokkos.
//
// data:  interleaved real/imag array of 2*n elements (1-based DFMatrixKokkos).
// n:     number of complex elements (= half the real signal length).
// isign: +1 forward (real -> half-complex), -1 inverse.
//
// Parallelization:
//   Twiddle factors for the recombination loop are precomputed, allowing
//   a FOR_ALL over all n/2 independent conjugate pairs.

#pragma once
#include <cmath>
#include "four1.hpp"

inline void realft(DFMatrixKokkos<double>& data, int n, int isign)
{
    double theta = 6.28318530717959 / (2.0 * n);
    double c1 = 0.5;
    double c2;

    if (isign == 1) {
        c2 = -0.5;
        four1(data, n, 1);
    } else {
        c2 = 0.5;
        theta = -theta;
    }

    int n2p3   = 2 * n + 3;
    int half_n = n / 2;

    // Precompute twiddle factors: index k corresponds to loop variable I = k+2.
    // Twiddle for I is exp(i*(I-1)*theta), so angle = (k+1)*theta.
    CArrayKokkos<double> tw_r(half_n), tw_i(half_n);
    FOR_ALL(k, 0, half_n, {
        double angle = (k + 1) * theta;
        tw_r(k) = cos(angle);
        tw_i(k) = sin(angle);
    });
    MATAR_FENCE();

    // Parallel recombination: each k touches indices (i1,i2) and (i3,i4)
    // which are in non-overlapping ranges across threads.
    FOR_ALL(k, 0, half_n, {
        int i  = k + 2;
        int i1 = 2 * i - 1;
        int i2 = i1 + 1;
        int i3 = n2p3 - i2;
        int i4 = i3 + 1;
        double wrs = tw_r(k);
        double wis = tw_i(k);
        double h1r =  c1 * (data(i1) + data(i3));
        double h1i =  c1 * (data(i2) - data(i4));
        double h2r = -c2 * (data(i2) + data(i4));
        double h2i =  c2 * (data(i1) - data(i3));
        data(i1) =  h1r + wrs * h2r - wis * h2i;
        data(i2) =  h1i + wrs * h2i + wis * h2r;
        data(i3) =  h1r - wrs * h2r + wis * h2i;
        data(i4) = -h1i + wrs * h2i + wis * h2r;
    });
    MATAR_FENCE();

    // DC / Nyquist fixup (single-element device operation)
    RUN({
        if (isign == 1) {
            double h1r = data(1);
            data(1) = h1r + data(2);
            data(2) = h1r - data(2);
        } else {
            double h1r = data(1);
            data(1) = c1 * (h1r + data(2));
            data(2) = c1 * (h1r - data(2));
        }
    });
    MATAR_FENCE();

    if (isign != 1) {
        four1(data, n, -1);
    }
}
