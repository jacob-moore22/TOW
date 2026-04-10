// Cubic spline evaluation (Numerical Recipes SPLINT).
// Given arrays xa[1..n], ya[1..n], and second derivatives y2a[1..n]
// (from spline()), evaluates the cubic spline at x, returning y.

#pragma once
#include <cstdio>
#include <matar.h>

using namespace mtr;

inline void splint(DFMatrixKokkos<double>& xa, DFMatrixKokkos<double>& ya,
                   DFMatrixKokkos<double>& y2a, int n, double x, double& y)
{
    int klo = 1;
    int khi = n;
    while (khi - klo > 1) {
        int k = (khi + klo) / 2;
        if (xa.host(k) > x) {
            khi = k;
        } else {
            klo = k;
        }
    }
    double h = xa.host(khi) - xa.host(klo);
    if (h == 0.0) {
        std::fprintf(stderr, "splint: bad xa input (identical values)\n");
        return;
    }
    double a = (xa.host(khi) - x) / h;
    double b = (x - xa.host(klo)) / h;
    y = a * ya.host(klo) + b * ya.host(khi) +
        ((a * a * a - a) * y2a.host(klo) +
         (b * b * b - b) * y2a.host(khi)) * (h * h) / 6.0;
}
