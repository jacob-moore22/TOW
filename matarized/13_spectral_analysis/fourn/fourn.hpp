// Parallel multi-dimensional FFT using MATAR + Kokkos.
//
// data:    interleaved real/imag array of 2*NTOT elements (1-based DFMatrixKokkos).
// nn_dims: host array of dimension sizes (0-based, nn_dims[0..ndim-1]).
//          Each dimension must be a power of 2.
// ndim:    number of dimensions.
// isign:   +1 forward, -1 inverse.
//
// Parallelization:
//   For each dimension: parallel bit-reversal (direct index computation)
//   and parallel butterfly stages with precomputed twiddle factors.

#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

inline void fourn(DFMatrixKokkos<double>& data, const int* nn_dims,
                  int ndim, int isign)
{
    int ntot = 1;
    for (int idim = 0; idim < ndim; idim++) ntot *= nn_dims[idim];

    int nprev = 1;
    for (int idim = 0; idim < ndim; idim++) {
        int n_dim = nn_dims[idim];
        int nrem  = ntot / (n_dim * nprev);
        int ip1   = 2 * nprev;
        int ip2   = ip1 * n_dim;

        // ---- Parallel bit-reversal along this dimension ----
        int log2_n = 0;
        for (int tmp = n_dim; tmp > 1; tmp >>= 1) log2_n++;

        int n_inner = nprev;
        int n_outer = nrem;
        int br_total = n_dim * n_inner * n_outer;

        int ip1_v = ip1, ip2_v = ip2, log2_v = log2_n;
        int ni = n_inner, no = n_outer;

        FOR_ALL(flat, 0, br_total, {
            int ci    = flat / (ni * no);
            int rem   = flat % (ni * no);
            int inner = rem / no;
            int outer = rem % no;

            int rev = 0;
            int tmp = ci;
            for (int b = 0; b < log2_v; b++) {
                rev = (rev << 1) | (tmp & 1);
                tmp >>= 1;
            }
            if (rev > ci) {
                int i = 1 + ci  * ip1_v + 2 * inner + outer * ip2_v;
                int j = 1 + rev * ip1_v + 2 * inner + outer * ip2_v;
                double tr = data(i);
                double ti = data(i + 1);
                data(i)     = data(j);
                data(i + 1) = data(j + 1);
                data(j)     = tr;
                data(j + 1) = ti;
            }
        });
        MATAR_FENCE();

        // ---- Butterfly stages along this dimension ----
        int max_tw = (n_dim > 1) ? n_dim / 2 : 1;
        CArrayKokkos<double> tw_r(max_tw), tw_i(max_tw);

        int ifp1 = ip1;
        while (ifp1 < ip2) {
            int ifp2 = 2 * ifp1;
            int half = ifp1 / ip1;
            double theta_base = isign * 6.28318530717959 / (ifp2 / ip1);

            FOR_ALL(mm, 0, half, {
                double angle = mm * theta_base;
                tw_r(mm) = cos(angle);
                tw_i(mm) = sin(angle);
            });
            MATAR_FENCE();

            int ngroups = 2 * ntot / ifp2;
            int total   = half * n_inner * ngroups;
            int ifp1_v  = ifp1, ifp2_v = ifp2;

            FOR_ALL(flat, 0, total, {
                int mm    = flat / (ni * ngroups);
                int rem2  = flat % (ni * ngroups);
                int inner2 = rem2 / ngroups;
                int group  = rem2 % ngroups;

                int k1 = 1 + mm * ip1_v + 2 * inner2 + group * ifp2_v;
                int k2 = k1 + ifp1_v;

                double wr = tw_r(mm);
                double wi = tw_i(mm);
                double tempr = wr * data(k2)     - wi * data(k2 + 1);
                double tempi = wr * data(k2 + 1) + wi * data(k2);
                data(k2)     = data(k1)     - tempr;
                data(k2 + 1) = data(k1 + 1) - tempi;
                data(k1)     = data(k1)     + tempr;
                data(k1 + 1) = data(k1 + 1) + tempi;
            });
            MATAR_FENCE();

            ifp1 = ifp2;
        }

        nprev *= n_dim;
    }
}
