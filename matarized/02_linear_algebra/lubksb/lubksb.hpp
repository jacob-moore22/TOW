// LU back-substitution (Numerical Recipes LUBKSB).
//
// Solves A·x = b given the LU decomposition of A (from ludcmp).
// a(1..n,1..n) is the LU matrix from ludcmp, indx(1..n) is the pivot
// vector, and b(1..n) is the RHS on input / solution on output.

#pragma once
#include <matar.h>

using namespace mtr;

inline void lubksb(DFMatrixKokkos<double>& a, int n,
                   DFMatrixKokkos<int>& indx, DFMatrixKokkos<double>& b)
{
    int ii = 0;
    for (int i = 1; i <= n; i++) {
        int ll = indx.host(i);
        double sum = b.host(ll);
        b.host(ll) = b.host(i);
        if (ii != 0) {
            for (int j = ii; j <= i - 1; j++)
                sum -= a.host(i, j) * b.host(j);
        } else if (sum != 0.0) {
            ii = i;
        }
        b.host(i) = sum;
    }

    for (int i = n; i >= 1; i--) {
        double sum = b.host(i);
        for (int j = i + 1; j <= n; j++)
            sum -= a.host(i, j) * b.host(j);
        b.host(i) = sum / a.host(i, i);
    }
}
