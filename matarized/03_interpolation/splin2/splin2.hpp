// 2D spline evaluation (Numerical Recipes SPLIN2).
// Given grids x1a[1..m], x2a[1..n], tabulated values ya(m,n), and
// second-derivative array y2a(m,n) (from splie2), evaluates the
// bicubic spline at (x1, x2) returning y.
// Depends on: spline, splint.

#pragma once
#include <matar.h>
#include "spline.hpp"
#include "splint.hpp"

using namespace mtr;

inline void splin2(DFMatrixKokkos<double>& x1a, DFMatrixKokkos<double>& x2a,
                   DFMatrixKokkos<double>& ya, DFMatrixKokkos<double>& y2a,
                   int m, int n,
                   double x1, double x2, double& y)
{
    DFMatrixKokkos<double> ytmp(n);
    DFMatrixKokkos<double> y2tmp(n);
    DFMatrixKokkos<double> yytmp(m);

    for (int j = 1; j <= m; j++) {
        for (int k = 1; k <= n; k++) {
            ytmp.host(k)  = ya.host(j, k);
            y2tmp.host(k) = y2a.host(j, k);
        }
        double tmp_y;
        splint(x2a, ytmp, y2tmp, n, x2, tmp_y);
        yytmp.host(j) = tmp_y;
    }

    DFMatrixKokkos<double> yy2tmp(m);
    spline(x1a, yytmp, m, 1.0e30, 1.0e30, yy2tmp);
    splint(x1a, yytmp, yy2tmp, m, x1, y);
}
