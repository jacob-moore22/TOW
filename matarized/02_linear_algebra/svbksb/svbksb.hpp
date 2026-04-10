// SVD back-substitution (Numerical Recipes SVBKSB).
//
// Solves A·x = b for x, given the SVD of A:  U(m,n), W(n), V(n,n).
// b(1..m) is the RHS, x(1..n) is the output solution.
// Small singular values should be zeroed in w before calling.

#pragma once
#include <matar.h>

using namespace mtr;

inline void svbksb(DFMatrixKokkos<double>& u, DFMatrixKokkos<double>& w,
                   DFMatrixKokkos<double>& v, int m, int n,
                   DFMatrixKokkos<double>& b, DFMatrixKokkos<double>& x)
{
    double tmp[n];

    for (int j = 1; j <= n; j++) {
        double s = 0.0;
        if (w.host(j) != 0.0) {
            for (int i = 1; i <= m; i++)
                s += u.host(i, j) * b.host(i);
            s /= w.host(j);
        }
        tmp[j - 1] = s;
    }

    for (int j = 1; j <= n; j++) {
        double s = 0.0;
        for (int jj = 1; jj <= n; jj++)
            s += v.host(j, jj) * tmp[jj - 1];
        x.host(j) = s;
    }
}
