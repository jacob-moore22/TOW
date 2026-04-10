// Parallel correlation of two real data sets via FFT using MATAR + Kokkos.
//
// data1, data2: real input arrays of n elements (1-based DFMatrixKokkos).
// n:            number of elements (must be a power of 2).
// ans:          output array of 2*n elements; after return ans(1..n)
//               holds the real-valued correlation.
//
// Parallelization:
//   Packing/unpacking handled by twofft.  Pointwise conjugate multiply
//   dispatched via DO_ALL over n/2+1 independent complex bins.

#pragma once
#include "twofft.hpp"
#include "realft.hpp"

inline void correl(DFMatrixKokkos<double>& data1,
                   DFMatrixKokkos<double>& data2,
                   int n,
                   DFMatrixKokkos<double>& ans)
{
    DFMatrixKokkos<double> fft(2 * n);

    twofft(data1, data2, fft, ans, n);

    int no2 = n / 2;

    // Pointwise: ANS(i) = FFT(i) * conj(ANS(i)) / no2
    DO_ALL(i, 1, no2 + 1, {
        double fr = fft(2 * i - 1);
        double fi = fft(2 * i);
        double ar = ans(2 * i - 1);
        double ai = ans(2 * i);
        ans(2 * i - 1) = (fr * ar + fi * ai) / no2;
        ans(2 * i)     = (fi * ar - fr * ai) / no2;
    });
    MATAR_FENCE();

    // Pack DC / Nyquist for realft inverse format
    RUN({
        ans(2) = ans(2 * no2 + 1);
    });
    MATAR_FENCE();

    realft(ans, no2, -1);
}
