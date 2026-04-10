// Tridiagonal system solver (Numerical Recipes TRIDAG).
//
// Solves a tridiagonal system of n equations:
//   a(i)*u(i-1) + b(i)*u(i) + c(i)*u(i+1) = r(i)
// where a(1) and c(n) are unused.
//
// All arrays are 1-based DFMatrixKokkos<double> of length n.
// u(1..n) contains the solution on return.

#pragma once
#include <cstdio>
#include <matar.h>

using namespace mtr;

inline void tridag(DFMatrixKokkos<double>& a, DFMatrixKokkos<double>& b,
                   DFMatrixKokkos<double>& c, DFMatrixKokkos<double>& r,
                   DFMatrixKokkos<double>& u, int n)
{
    double gam[n];

    if (b.host(1) == 0.0) {
        std::fprintf(stderr, "tridag: zero pivot at row 1\n");
        return;
    }

    double bet = b.host(1);
    u.host(1) = r.host(1) / bet;

    for (int j = 2; j <= n; j++) {
        gam[j - 1] = c.host(j - 1) / bet;
        bet = b.host(j) - a.host(j) * gam[j - 1];
        if (bet == 0.0) {
            std::fprintf(stderr, "tridag: zero pivot at row %d\n", j);
            return;
        }
        u.host(j) = (r.host(j) - a.host(j) * u.host(j - 1)) / bet;
    }

    for (int j = n - 1; j >= 1; j--)
        u.host(j) -= gam[j] * u.host(j + 1);
}
