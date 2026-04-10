#pragma once
#include <algorithm>
#include <matar.h>

using namespace mtr;

// Evaluate polynomial c(1..nc) and all its derivatives at point x.
// On output, pd(i) = (i-1)th derivative of the polynomial at x.
// Coefficients: c(1) = constant, c(2) = x^1, ..., c(nc) = x^(nc-1).
inline void ddpoly(const DFMatrixKokkos<double>& c, int nc, double x,
                   DFMatrixKokkos<double>& pd, int nd)
{
    pd(1) = c(nc);
    for (int j = 2; j <= nd; j++)
        pd(j) = 0.0;

    for (int i = nc - 1; i >= 1; i--) {
        int nnd = std::min(nd, nc + 1 - i);
        for (int j = nnd; j >= 2; j--)
            pd(j) = pd(j) * x + pd(j - 1);
        pd(1) = pd(1) * x + c(i);
    }

    double cnst = 2.0;
    for (int i = 3; i <= nd; i++) {
        pd(i) *= cnst;
        cnst  *= i;
    }
}
