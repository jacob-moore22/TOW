// Cubic spline coefficient computation (Numerical Recipes SPLINE).
// Given arrays x[1..n] and y[1..n] of tabulated function values, and
// first derivatives yp1 and ypn at the endpoints (set > 0.99e30 for
// natural spline), computes second derivatives y2[1..n].

#pragma once
#include <matar.h>

using namespace mtr;

inline void spline(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                   int n, double yp1, double ypn,
                   DFMatrixKokkos<double>& y2)
{
    constexpr int NMAX = 100;
    double u[NMAX + 1];

    if (yp1 > 0.99e30) {
        y2.host(1) = 0.0;
        u[1]       = 0.0;
    } else {
        y2.host(1) = -0.5;
        u[1] = (3.0 / (x.host(2) - x.host(1))) *
               ((y.host(2) - y.host(1)) / (x.host(2) - x.host(1)) - yp1);
    }

    for (int i = 2; i <= n - 1; i++) {
        double sig = (x.host(i) - x.host(i - 1)) /
                     (x.host(i + 1) - x.host(i - 1));
        double p   = sig * y2.host(i - 1) + 2.0;
        y2.host(i) = (sig - 1.0) / p;
        u[i] = (6.0 * ((y.host(i + 1) - y.host(i)) /
                        (x.host(i + 1) - x.host(i)) -
                        (y.host(i) - y.host(i - 1)) /
                        (x.host(i) - x.host(i - 1))) /
                (x.host(i + 1) - x.host(i - 1)) - sig * u[i - 1]) / p;
    }

    double qn, un;
    if (ypn > 0.99e30) {
        qn = 0.0;
        un = 0.0;
    } else {
        qn = 0.5;
        un = (3.0 / (x.host(n) - x.host(n - 1))) *
             (ypn - (y.host(n) - y.host(n - 1)) /
              (x.host(n) - x.host(n - 1)));
    }

    y2.host(n) = (un - qn * u[n - 1]) / (qn * y2.host(n - 1) + 1.0);

    for (int k = n - 1; k >= 1; k--) {
        y2.host(k) = y2.host(k) * y2.host(k + 1) + u[k];
    }
}
