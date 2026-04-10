// Rational function interpolation (Numerical Recipes RATINT).
// Given arrays xa[1..n] and ya[1..n], and a value x, returns y and
// error estimate dy using the diagonal rational function method.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

inline void ratint(DFMatrixKokkos<double>& xa, DFMatrixKokkos<double>& ya,
                   int n, double x, double& y, double& dy)
{
    constexpr int    NMAX = 10;
    constexpr double TINY = 1.0e-25;
    double c[NMAX + 1], d[NMAX + 1];

    int ns = 1;
    double hh = std::fabs(x - xa.host(1));
    for (int i = 1; i <= n; i++) {
        double h = std::fabs(x - xa.host(i));
        if (h == 0.0) {
            y  = ya.host(i);
            dy = 0.0;
            return;
        } else if (h < hh) {
            ns = i;
            hh = h;
        }
        c[i] = ya.host(i);
        d[i] = ya.host(i) + TINY;
    }
    y  = ya.host(ns);
    ns = ns - 1;
    for (int m = 1; m <= n - 1; m++) {
        for (int i = 1; i <= n - m; i++) {
            double w  = c[i + 1] - d[i];
            double h  = xa.host(i + m) - x;
            double t  = (xa.host(i) - x) * d[i] / h;
            double dd = t - c[i + 1];
            if (dd == 0.0) {
                std::fprintf(stderr, "ratint: pole encountered\n");
                return;
            }
            dd   = w / dd;
            d[i] = c[i + 1] * dd;
            c[i] = t * dd;
        }
        if (2 * ns < n - m) {
            dy = c[ns + 1];
        } else {
            dy = d[ns];
            ns = ns - 1;
        }
        y = y + dy;
    }
}
