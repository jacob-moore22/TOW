// Parallel convolution / deconvolution via FFT using MATAR + Kokkos.
//
// data:   real signal array of n elements (1-based DFMatrixKokkos, device).
// n:      signal length (must be power of 2).
// respns: real response array of n elements (modified: padded in-place).
// m:      original response length (odd, < n). Response is stored in
//         respns(1..m); the rest is used as scratch for wrap-around padding.
// isign:  +1 convolution, -1 deconvolution.
// ans:    output array of 2*n elements. After return, ans(1..n) contains
//         the real-valued convolution result.
//
// Parallelization:
//   Wrap-around copy  -- DO_ALL over (m-1)/2 elements.
//   Zero fill         -- DO_ALL over the middle gap.
//   Pointwise multiply -- DO_ALL over n/2+1 complex bins (embarrassingly parallel).
//   FFT stages are parallelized internally by four1 / realft / twofft.

#pragma once
#include "twofft.hpp"
#include "realft.hpp"

inline void convlv(DFMatrixKokkos<double>& data,
                   int n,
                   DFMatrixKokkos<double>& respns,
                   int m,
                   int isign,
                   DFMatrixKokkos<double>& ans)
{
    // Wrap-around: copy response tail to end of padded array
    int half_m = (m - 1) / 2;
    if (half_m > 0) {
        DO_ALL(i, 1, half_m, {
            respns(n + 1 - i) = respns(m + 1 - i);
        });
        MATAR_FENCE();
    }

    // Zero-fill the middle of the padded response
    int zero_start = (m + 3) / 2;
    int zero_end   = n - half_m;
    if (zero_start <= zero_end) {
        DO_ALL(i, zero_start, zero_end, {
            respns(i) = 0.0;
        });
        MATAR_FENCE();
    }

    // FFT work array (n complex elements = 2n reals)
    DFMatrixKokkos<double> fft(2 * n);

    twofft(data, respns, fft, ans, n);

    int no2 = n / 2;

    // Pointwise complex multiply (convolution) or divide (deconvolution)
    DO_ALL(i, 1, no2 + 1, {
        double fr = fft(2 * i - 1);
        double fi = fft(2 * i);
        double ar = ans(2 * i - 1);
        double ai = ans(2 * i);
        if (isign == 1) {
            ans(2 * i - 1) = (fr * ar - fi * ai) / no2;
            ans(2 * i)     = (fr * ai + fi * ar) / no2;
        } else {
            double denom = (ar * ar + ai * ai) * no2;
            ans(2 * i - 1) = (fr * ar + fi * ai) / denom;
            ans(2 * i)     = (fi * ar - fr * ai) / denom;
        }
    });
    MATAR_FENCE();

    // Pack DC and Nyquist components for REALFT inverse format:
    // ans(2) ← real part of ans(no2+1)
    RUN({
        ans(2) = ans(2 * no2 + 1);
    });
    MATAR_FENCE();

    realft(ans, no2, -1);
}
