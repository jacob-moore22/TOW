// Iterative improvement of a solution vector (Numerical Recipes MPROVE).
//
// a(1..n,1..n)    -- original coefficient matrix (unmodified)
// alud(1..n,1..n) -- LU decomposition of a (from ludcmp)
// indx(1..n)      -- pivot indices (from ludcmp)
// b(1..n)         -- original RHS vector
// x(1..n)         -- solution vector (improved in-place)
//
// Uses double-precision accumulation for the residual to recover
// accuracy lost in the original solve.

#pragma once
#include <cmath>
#include <matar.h>
#include "lubksb.hpp"

using namespace mtr;

inline void mprove(DFMatrixKokkos<double>& a, DFMatrixKokkos<double>& alud,
                   int n, DFMatrixKokkos<int>& indx,
                   DFMatrixKokkos<double>& b, DFMatrixKokkos<double>& x)
{
    DFMatrixKokkos<double> r(n);

    for (int i = 1; i <= n; i++) {
        double sdp = -static_cast<double>(b.host(i));
        for (int j = 1; j <= n; j++)
            sdp += static_cast<double>(a.host(i, j)) * static_cast<double>(x.host(j));
        r.host(i) = sdp;
    }

    lubksb(alud, n, indx, r);

    for (int i = 1; i <= n; i++)
        x.host(i) -= r.host(i);
}
