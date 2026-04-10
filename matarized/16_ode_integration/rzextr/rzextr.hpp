// Rational function extrapolation for Bulirsch-Stoer method
// (Numerical Recipes RZEXTR).
// Called iteratively with iest=1,2,... to build up an extrapolation tableau.
// Static arrays x[] and d[][] persist between calls within a BS step sequence.

#pragma once
#include <cmath>
#include <algorithm>
#include <matar.h>

using namespace mtr;

inline void rzextr(int iest, double xest, DFMatrixKokkos<double>& yest,
                   DFMatrixKokkos<double>& yz, DFMatrixKokkos<double>& dy,
                   int nv, int nuse)
{
    constexpr int IMAX = 11;
    constexpr int NMAX = 10;
    constexpr int NCOL = 7;

    static double x[IMAX + 1];
    static double d[NMAX + 1][NCOL + 1];
    double fx[NCOL + 1];

    x[iest] = xest;

    if (iest == 1) {
        for (int j = 1; j <= nv; j++) {
            yz.host(j) = yest.host(j);
            d[j][1]    = yest.host(j);
            dy.host(j) = yest.host(j);
        }
    } else {
        int m1 = std::min(iest, nuse);
        for (int k = 1; k <= m1 - 1; k++)
            fx[k + 1] = x[iest - k] / xest;

        for (int j = 1; j <= nv; j++) {
            double yy  = yest.host(j);
            double v   = d[j][1];
            double c   = yy;
            d[j][1]    = yy;
            double ddy = 0.0;

            for (int k = 2; k <= m1; k++) {
                double b1 = fx[k] * v;
                double b  = b1 - c;
                if (b != 0.0) {
                    b   = (c - v) / b;
                    ddy = c * b;
                    c   = b1 * b;
                } else {
                    ddy = v;
                }
                v      = d[j][k];
                d[j][k] = ddy;
                yy    += ddy;
            }
            dy.host(j) = ddy;
            yz.host(j) = yy;
        }
    }
}
