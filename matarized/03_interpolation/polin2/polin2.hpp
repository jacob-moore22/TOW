// 2D polynomial interpolation (Numerical Recipes POLIN2).
// Given grids x1a[1..m], x2a[1..n], and tabulated values ya(m,n),
// returns interpolated value y and error estimate dy at (x1, x2).
// Uses polint for row-wise and then column-wise interpolation.

#pragma once
#include <matar.h>
#include "polint.hpp"

using namespace mtr;

inline void polin2(DFMatrixKokkos<double>& x1a, DFMatrixKokkos<double>& x2a,
                   DFMatrixKokkos<double>& ya, int m, int n,
                   double x1, double x2, double& y, double& dy)
{
    constexpr int NMAX = 20;
    DFMatrixKokkos<double> yntmp(NMAX);
    DFMatrixKokkos<double> ymtmp(NMAX);

    for (int j = 1; j <= m; j++) {
        for (int k = 1; k <= n; k++) {
            yntmp.host(k) = ya.host(j, k);
        }
        double tmp_y, tmp_dy;
        polint(x2a, yntmp, n, x2, tmp_y, tmp_dy);
        ymtmp.host(j) = tmp_y;
    }
    polint(x1a, ymtmp, m, x1, y, dy);
}
