// 2D spline setup (Numerical Recipes SPLIE2).
// Given grids x1a[1..m], x2a[1..n], and tabulated values ya(m,n),
// computes second-derivative array y2a(m,n) by applying 1D natural
// cubic spline along each row (x2 direction).
// Depends on: spline.

#pragma once
#include <matar.h>
#include "spline.hpp"

using namespace mtr;

inline void splie2(DFMatrixKokkos<double>& x1a, DFMatrixKokkos<double>& x2a,
                   DFMatrixKokkos<double>& ya, int m, int n,
                   DFMatrixKokkos<double>& y2a)
{
    DFMatrixKokkos<double> ytmp(n);
    DFMatrixKokkos<double> y2tmp(n);

    for (int j = 1; j <= m; j++) {
        for (int k = 1; k <= n; k++) {
            ytmp.host(k) = ya.host(j, k);
        }
        spline(x2a, ytmp, n, 1.0e30, 1.0e30, y2tmp);
        for (int k = 1; k <= n; k++) {
            y2a.host(j, k) = y2tmp.host(k);
        }
    }
}
