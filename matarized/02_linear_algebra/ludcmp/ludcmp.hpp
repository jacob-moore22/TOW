// LU decomposition with partial pivoting (Numerical Recipes LUDCMP).
//
// Given a(1..n, 1..n), replaces it with the LU decomposition of a
// rowwise permutation of itself. indx(1..n) records the row permutation.
// d is +/-1 depending on whether the number of row interchanges was
// even or odd.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

inline void ludcmp(DFMatrixKokkos<double>& a, int n,
                   DFMatrixKokkos<int>& indx, double& d)
{
    constexpr double TINY = 1.0e-20;
    double vv[n];

    d = 1.0;
    for (int i = 1; i <= n; i++) {
        double aamax = 0.0;
        for (int j = 1; j <= n; j++) {
            double tmp = std::fabs(a.host(i, j));
            if (tmp > aamax) aamax = tmp;
        }
        if (aamax == 0.0) {
            std::fprintf(stderr, "ludcmp: singular matrix\n");
            return;
        }
        vv[i - 1] = 1.0 / aamax;
    }

    for (int j = 1; j <= n; j++) {
        for (int i = 1; i <= j - 1; i++) {
            double sum = a.host(i, j);
            for (int k = 1; k <= i - 1; k++)
                sum -= a.host(i, k) * a.host(k, j);
            a.host(i, j) = sum;
        }

        double aamax = 0.0;
        int imax = j;
        for (int i = j; i <= n; i++) {
            double sum = a.host(i, j);
            for (int k = 1; k <= j - 1; k++)
                sum -= a.host(i, k) * a.host(k, j);
            a.host(i, j) = sum;
            double dum = vv[i - 1] * std::fabs(sum);
            if (dum >= aamax) {
                imax = i;
                aamax = dum;
            }
        }

        if (j != imax) {
            for (int k = 1; k <= n; k++) {
                double dum = a.host(imax, k);
                a.host(imax, k) = a.host(j, k);
                a.host(j, k) = dum;
            }
            d = -d;
            vv[imax - 1] = vv[j - 1];
        }

        indx.host(j) = imax;

        if (j != n) {
            if (a.host(j, j) == 0.0) a.host(j, j) = TINY;
            double dum = 1.0 / a.host(j, j);
            for (int i = j + 1; i <= n; i++)
                a.host(i, j) *= dum;
        }
    }
    if (a.host(n, n) == 0.0) a.host(n, n) = TINY;
}
