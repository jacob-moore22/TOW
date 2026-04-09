// Parallel simultaneous FFT of two real sequences using MATAR + Kokkos.
//
// data1, data2: real input arrays of n elements (1-based DFMatrixKokkos).
// fft1, fft2:   interleaved complex output arrays of 2*n elements each.
//               Complex element k: real at fft(2k-1), imag at fft(2k).
// n:            number of real elements (must be power of 2).
//
// Parallelization:
//   Pack:   embarrassingly parallel over n elements.
//   Unpack: conjugate-symmetry separation -- each thread j writes to
//           positions j and n+2-j, which are non-overlapping across threads.

#pragma once
#include "four1.hpp"

inline void twofft(DFMatrixKokkos<double>& data1,
                   DFMatrixKokkos<double>& data2,
                   DFMatrixKokkos<double>& fft1,
                   DFMatrixKokkos<double>& fft2,
                   int n)
{
    // Pack two real sequences into one complex sequence
    DO_ALL(j, 1, n, {
        fft1(2 * j - 1) = data1(j);
        fft1(2 * j)     = data2(j);
    });
    MATAR_FENCE();

    four1(fft1, n, 1);

    // Separate the DC components
    RUN({
        double re1 = fft1(1);
        double im1 = fft1(2);
        fft2(1) = im1;
        fft2(2) = 0.0;
        fft1(1) = re1;
        fft1(2) = 0.0;
    });
    MATAR_FENCE();

    int n2 = n + 2;

    // Unpack: separate two spectra via conjugate symmetry.
    // Thread j accesses positions j and n2-j; these pairs are disjoint
    // across different j values (j in [2, n/2+1], n2-j in [n/2+1, n]).
    DO_ALL(j, 2, n / 2 + 1, {
        int kk = n2 - j;
        double f_jr = fft1(2 * j - 1);
        double f_ji = fft1(2 * j);
        double f_kr = fft1(2 * kk - 1);
        double f_ki = fft1(2 * kk);

        // H1 = (0.5, 0) * (FFT1(j) + conj(FFT1(kk)))
        double h1r = 0.5 * (f_jr + f_kr);
        double h1i = 0.5 * (f_ji - f_ki);

        // H2 = (0, -0.5) * (FFT1(j) - conj(FFT1(kk)))
        double h2r =  0.5 * (f_ji + f_ki);
        double h2i = -0.5 * (f_jr - f_kr);

        fft1(2 * j - 1)  =  h1r;
        fft1(2 * j)      =  h1i;
        fft1(2 * kk - 1) =  h1r;
        fft1(2 * kk)     = -h1i;
        fft2(2 * j - 1)  =  h2r;
        fft2(2 * j)      =  h2i;
        fft2(2 * kk - 1) =  h2r;
        fft2(2 * kk)     = -h2i;
    });
    MATAR_FENCE();
}
