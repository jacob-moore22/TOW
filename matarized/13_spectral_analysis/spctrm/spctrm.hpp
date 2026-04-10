// Parallel power-spectrum estimation (Bartlett window) using MATAR + Kokkos.
//
// p:          output power spectrum, m elements (1-based DFMatrixKokkos).
// m:          number of spectral estimates.
// k:          number of data segments to average.
// ovrlap:     true to overlap successive segments.
// w1:         work array of at least 4*m elements (device).
// w2:         work array of m elements (device).
// input_data: flat data array (0-based DCArrayKokkos, all segments sequentially).
// read_pos:   in/out index into input_data (0-based, advanced by the routine).
//
// Parallelization:
//   Windowing -- DO_ALL over 2*m independent pairs.
//   Power accumulation -- DO_ALL over m-1 independent bins.

#pragma once
#include <cmath>
#include <matar.h>
#include "four1.hpp"

using namespace mtr;

inline void spctrm(DFMatrixKokkos<double>& p, int m, int k, bool ovrlap,
                   DFMatrixKokkos<double>& w1,
                   DFMatrixKokkos<double>& w2,
                   DCArrayKokkos<double>& input_data, int& read_pos)
{
    int mm  = 2 * m;
    int m4  = 4 * m;
    int m44 = m4 + 4;
    int m43 = m4 + 3;

    double facm = m - 0.5;
    double facp = 1.0 / (m + 0.5);

    // Sum of squared Bartlett window weights
    double sumw = 0.0;
    double loc_sw = 0.0;
    FOR_REDUCE_SUM(j, 0, mm, loc_sw, {
        double wval = 1.0 - fabs((j - facm) * facp);
        loc_sw += wval * wval;
    }, sumw);

    DO_ALL(j, 1, m, { p(j) = 0.0; });
    MATAR_FENCE();

    // Helper: read m values from input_data into a DFMatrixKokkos (host→device)
    auto read_m = [&](DFMatrixKokkos<double>& buf) {
        for (int j = 1; j <= m; j++)
            buf.host(j) = input_data.host(read_pos++);
        buf.update_device();
    };

    if (ovrlap) read_m(w2);

    double den = 0.0;
    for (int kk = 1; kk <= k; kk++) {
        if (ovrlap) {
            for (int joff = -1; joff <= 0; joff++) {
                int joff_v = joff;
                DO_ALL(j, 1, m, { w1(joff_v + 2 * j) = w2(j); });
                MATAR_FENCE();

                read_m(w2);

                int joffn = joff + mm;
                DO_ALL(j, 1, m, { w1(joffn + 2 * j) = w2(j); });
                MATAR_FENCE();
            }
        } else {
            for (int joff = -1; joff <= 0; joff++) {
                int start = joff + 2;
                for (int pos = start; pos <= m4; pos += 2)
                    w1.host(pos) = input_data.host(read_pos++);
            }
            w1.update_device();
        }

        // Apply Bartlett window (parallel)
        DO_ALL(j, 1, mm, {
            int j2 = 2 * j;
            double w = 1.0 - fabs(((j - 1) - facm) * facp);
            w1(j2)     *= w;
            w1(j2 - 1) *= w;
        });
        MATAR_FENCE();

        four1(w1, mm, 1);

        // Accumulate power spectrum
        RUN({ p(1) += w1(1) * w1(1) + w1(2) * w1(2); });
        MATAR_FENCE();

        DO_ALL(j, 2, m, {
            int j2 = 2 * j;
            p(j) += w1(j2) * w1(j2) + w1(j2 - 1) * w1(j2 - 1)
                  + w1(m44 - j2) * w1(m44 - j2)
                  + w1(m43 - j2) * w1(m43 - j2);
        });
        MATAR_FENCE();

        den += sumw;
    }

    double den_inv = 1.0 / (m4 * den);
    DO_ALL(j, 1, m, { p(j) *= den_inv; });
    MATAR_FENCE();
}
