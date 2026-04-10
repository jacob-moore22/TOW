// Polynomial extrapolation for Bulirsch-Stoer method (Numerical Recipes PZEXTR).
// Called iteratively with iest=1,2,... to build up an extrapolation tableau.
// Static arrays x[] and qcol[][] persist between calls within a BS step sequence.

#pragma once
#include <cmath>
#include <algorithm>
#include <matar.h>

using namespace mtr;

inline void pzextr(int iest, double xest, DFMatrixKokkos<double>& yest,
                   DFMatrixKokkos<double>& yz, DFMatrixKokkos<double>& dy,
                   int nv, int nuse)
{
    constexpr int IMAX = 11;
    constexpr int NCOL = 7;
    constexpr int NMAX = 10;

    static double x[IMAX + 1];
    static double qcol[NMAX + 1][NCOL + 1];
    double d[NMAX + 1];

    x[iest] = xest;
    for (int j = 1; j <= nv; j++) {
        dy.host(j) = yest.host(j);
        yz.host(j) = yest.host(j);
    }

    if (iest == 1) {
        for (int j = 1; j <= nv; j++)
            qcol[j][1] = yest.host(j);
    } else {
        int m1 = std::min(iest, nuse);
        for (int j = 1; j <= nv; j++)
            d[j] = yest.host(j);

        for (int k1 = 1; k1 <= m1 - 1; k1++) {
            double delta = 1.0 / (x[iest - k1] - xest);
            double f1 = xest * delta;
            double f2 = x[iest - k1] * delta;
            for (int j = 1; j <= nv; j++) {
                double q  = qcol[j][k1];
                qcol[j][k1] = dy.host(j);
                delta     = d[j] - q;
                dy.host(j) = f1 * delta;
                d[j]      = f2 * delta;
                yz.host(j) += dy.host(j);
            }
        }

        for (int j = 1; j <= nv; j++)
            qcol[j][m1] = dy.host(j);
    }
}
