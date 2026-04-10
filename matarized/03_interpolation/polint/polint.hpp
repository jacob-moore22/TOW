// Neville's polynomial interpolation (Numerical Recipes POLINT).
// Given arrays xa[1..n] and ya[1..n], and a value x, returns y and
// error estimate dy via Neville's algorithm.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

inline void polint(DFMatrixKokkos<double>& xa, DFMatrixKokkos<double>& ya,
                   int n, double x, double& y, double& dy)
{
    constexpr int NMAX = 10;
    double c[NMAX + 1], d[NMAX + 1];

    int ns = 1;
    double dif = std::fabs(x - xa.host(1));
    for (int i = 1; i <= n; i++) {
        double dift = std::fabs(x - xa.host(i));
        if (dift < dif) {
            ns  = i;
            dif = dift;
        }
        c[i] = ya.host(i);
        d[i] = ya.host(i);
    }
    y  = ya.host(ns);
    ns = ns - 1;
    for (int m = 1; m <= n - 1; m++) {
        for (int i = 1; i <= n - m; i++) {
            double ho  = xa.host(i) - x;
            double hp  = xa.host(i + m) - x;
            double w   = c[i + 1] - d[i];
            double den = ho - hp;
            if (den == 0.0) {
                std::fprintf(stderr, "polint: identical xa entries\n");
                return;
            }
            den  = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
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
