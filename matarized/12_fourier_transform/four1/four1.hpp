// Parallel Cooley-Tukey radix-2 FFT using MATAR + Kokkos.
//
// data: interleaved real/imag array of 2*nn elements (1-based DFMatrixKokkos).
//       Complex element k has real part at data(2k-1), imag at data(2k).
// nn:   number of complex elements (must be power of 2).
// isign: +1 forward, -1 inverse.
//
// Parallelization:
//   Phase 1 -- bit-reversal computed directly per index (no sequential recurrence).
//   Phase 2 -- butterfly stages: twiddle factors precomputed per stage,
//              then all butterflies dispatched via FOR_ALL.

#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

inline void four1(DFMatrixKokkos<double>& data, int nn, int isign)
{
    int n = 2 * nn;

    int log2_nn = 0;
    for (int tmp = nn; tmp > 1; tmp >>= 1) log2_nn++;

    // ---- Phase 1: parallel bit-reversal permutation ----
    // Each thread computes the bit-reversed partner of its complex index.
    // The guard rev > ci ensures each pair is swapped exactly once.
    FOR_ALL(ci, 0, nn, {
        int rev = 0;
        int tmp = ci;
        for (int b = 0; b < log2_nn; b++) {
            rev = (rev << 1) | (tmp & 1);
            tmp >>= 1;
        }
        if (rev > ci) {
            int i = 2 * ci + 1;
            int j = 2 * rev + 1;
            double tr = data(j);
            double ti = data(j + 1);
            data(j)     = data(i);
            data(j + 1) = data(i + 1);
            data(i)     = tr;
            data(i + 1) = ti;
        }
    });
    MATAR_FENCE();

    // ---- Phase 2: butterfly stages ----
    // Pre-allocate twiddle arrays at max size needed (nn/2).
    int max_tw = (nn > 1) ? nn / 2 : 1;
    CArrayKokkos<double> tw_r(max_tw), tw_i(max_tw);

    int mmax = 2;
    while (n > mmax) {
        int istep = 2 * mmax;
        int half  = mmax / 2;
        double theta_base = 6.28318530717959 / (isign * mmax);

        FOR_ALL(m, 0, half, {
            double angle = theta_base * m;
            tw_r(m) = cos(angle);
            tw_i(m) = sin(angle);
        });
        MATAR_FENCE();

        int ngroups = n / istep;
        int total   = ngroups * half;

        // Each k uniquely maps to one butterfly (group, twiddle index).
        // Within a stage, no two butterflies share array elements.
        FOR_ALL(k, 0, total, {
            int group = k / half;
            int m     = k % half;
            int i = 1 + group * istep + 2 * m;
            int j = i + mmax;
            double wr = tw_r(m);
            double wi = tw_i(m);
            double tempr = wr * data(j)     - wi * data(j + 1);
            double tempi = wr * data(j + 1) + wi * data(j);
            data(j)     = data(i)     - tempr;
            data(j + 1) = data(i + 1) - tempi;
            data(i)     = data(i)     + tempr;
            data(i + 1) = data(i + 1) + tempi;
        });
        MATAR_FENCE();

        mmax = istep;
    }
}
