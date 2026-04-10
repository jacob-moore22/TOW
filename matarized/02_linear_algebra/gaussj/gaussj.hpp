// Gauss-Jordan elimination with full pivoting (Numerical Recipes GAUSSJ).
//
// Solves A·X = B in-place. On return A is replaced by its inverse and
// B is replaced by the solution vectors.
//
// a(1..n, 1..n)  -- coefficient matrix (overwritten with inverse)
// n              -- logical dimension
// b(1..n, 1..m)  -- right-hand side matrix (overwritten with solutions)
// m              -- number of RHS columns

#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <matar.h>

using namespace mtr;

inline void gaussj(DFMatrixKokkos<double>& a, int n,
                   DFMatrixKokkos<double>& b, int m)
{
    int indxc[n], indxr[n], ipiv[n];

    for (int j = 1; j <= n; j++) ipiv[j - 1] = 0;

    for (int i = 1; i <= n; i++) {
        double big = 0.0;
        int irow = 0, icol = 0;

        for (int j = 1; j <= n; j++) {
            if (ipiv[j - 1] != 1) {
                for (int k = 1; k <= n; k++) {
                    if (ipiv[k - 1] == 0) {
                        if (std::fabs(a.host(j, k)) >= big) {
                            big  = std::fabs(a.host(j, k));
                            irow = j;
                            icol = k;
                        }
                    } else if (ipiv[k - 1] > 1) {
                        std::fprintf(stderr, "gaussj: singular matrix (pivot)\n");
                        return;
                    }
                }
            }
        }
        ipiv[icol - 1] += 1;

        if (irow != icol) {
            for (int l = 1; l <= n; l++) {
                double tmp = a.host(irow, l);
                a.host(irow, l) = a.host(icol, l);
                a.host(icol, l) = tmp;
            }
            for (int l = 1; l <= m; l++) {
                double tmp = b.host(irow, l);
                b.host(irow, l) = b.host(icol, l);
                b.host(icol, l) = tmp;
            }
        }

        indxr[i - 1] = irow;
        indxc[i - 1] = icol;

        if (a.host(icol, icol) == 0.0) {
            std::fprintf(stderr, "gaussj: singular matrix\n");
            return;
        }

        double pivinv = 1.0 / a.host(icol, icol);
        a.host(icol, icol) = 1.0;
        for (int l = 1; l <= n; l++) a.host(icol, l) *= pivinv;
        for (int l = 1; l <= m; l++) b.host(icol, l) *= pivinv;

        for (int ll = 1; ll <= n; ll++) {
            if (ll != icol) {
                double dum = a.host(ll, icol);
                a.host(ll, icol) = 0.0;
                for (int l = 1; l <= n; l++)
                    a.host(ll, l) -= a.host(icol, l) * dum;
                for (int l = 1; l <= m; l++)
                    b.host(ll, l) -= b.host(icol, l) * dum;
            }
        }
    }

    for (int l = n; l >= 1; l--) {
        if (indxr[l - 1] != indxc[l - 1]) {
            for (int k = 1; k <= n; k++) {
                double tmp = a.host(k, indxr[l - 1]);
                a.host(k, indxr[l - 1]) = a.host(k, indxc[l - 1]);
                a.host(k, indxc[l - 1]) = tmp;
            }
        }
    }
}
