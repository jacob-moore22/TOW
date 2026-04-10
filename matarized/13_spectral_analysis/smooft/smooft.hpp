// Parallel FFT-based data smoothing using MATAR + Kokkos.
//
// y:   real array of n elements (1-based DFMatrixKokkos, modified in-place).
// n:   signal length.
// pts: approximate smoothing width in data points (controls the cutoff
//      frequency of a low-pass Gaussian-like filter).
//
// Parallelization:
//   Trend removal / restoration -- DO_ALL over n elements.
//   Zero padding               -- DO_ALL.
//   Frequency-domain filter     -- DO_ALL over mo2-1 bins.

#pragma once
#include <cmath>
#include <matar.h>
#include "realft.hpp"

using namespace mtr;

inline void smooft(DFMatrixKokkos<double>& y, int n, double pts)
{
    // Find next power of 2 >= n + 2*pts
    int m = 2;
    int nmin = n + static_cast<int>(2.0 * pts);
    while (m < nmin) m *= 2;

    // Work array (m may exceed y's allocation)
    DFMatrixKokkos<double> w(m);

    // Capture y(1) and y(n) on device before modification
    CArrayKokkos<double> ep(2);
    RUN({ ep(0) = y(1); ep(1) = y(n); });
    MATAR_FENCE();

    double rn1 = 1.0 / (n - 1.0);

    // Remove linear trend and copy to work array
    DO_ALL(j, 1, n, {
        w(j) = y(j) - rn1 * (ep(0) * (n - j) + ep(1) * (j - 1));
    });
    MATAR_FENCE();

    // Zero-pad
    if (n + 1 <= m) {
        DO_ALL(j, n + 1, m, { w(j) = 0.0; });
        MATAR_FENCE();
    }

    int mo2 = m / 2;
    double cnst = (pts / m) * (pts / m);

    realft(w, mo2, 1);

    // DC component
    RUN({ w(1) = w(1) / mo2; });
    MATAR_FENCE();

    // Low-pass filter in frequency domain
    if (mo2 > 1) {
        DO_ALL(j, 1, mo2 - 1, {
            double raw = (1.0 - cnst * j * j) / mo2;
            double fac = (raw > 0.0) ? raw : 0.0;
            int k = 2 * j + 1;
            w(k)     *= fac;
            w(k + 1) *= fac;
        });
        MATAR_FENCE();
    }

    // Nyquist component
    double nyq_raw = (1.0 - 0.25 * pts * pts) / mo2;
    double nyq_fac = (nyq_raw > 0.0) ? nyq_raw : 0.0;
    RUN({ w(2) *= nyq_fac; });
    MATAR_FENCE();

    realft(w, mo2, -1);

    // Restore linear trend and write back
    DO_ALL(j, 1, n, {
        y(j) = rn1 * (ep(0) * (n - j) + ep(1) * (j - 1)) + w(j);
    });
    MATAR_FENCE();
}
