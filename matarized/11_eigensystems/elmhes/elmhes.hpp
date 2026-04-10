// Reduction to upper Hessenberg form by Gaussian elimination with pivoting
// (Numerical Recipes ELMHES).
// Replaces a[1..n][1..n] by an upper Hessenberg matrix with identical
// eigenvalues. Only elements a[i][j] with i <= j+1 are meaningful on output.
// Recommended preprocessing: balanc.

#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

inline void elmhes(DFMatrixKokkos<double>& a, int n)
{
    if (n > 2) {
        for (int m = 2; m <= n - 1; m++) {
            double x = 0.0;
            int i = m;
            for (int j = m; j <= n; j++) {
                if (std::fabs(a.host(j, m - 1)) > std::fabs(x)) {
                    x = a.host(j, m - 1);
                    i = j;
                }
            }
            if (i != m) {
                for (int j = m - 1; j <= n; j++) {
                    double y = a.host(i, j);
                    a.host(i, j) = a.host(m, j);
                    a.host(m, j) = y;
                }
                for (int j = 1; j <= n; j++) {
                    double y = a.host(j, i);
                    a.host(j, i) = a.host(j, m);
                    a.host(j, m) = y;
                }
            }
            if (x != 0.0) {
                for (i = m + 1; i <= n; i++) {
                    double y = a.host(i, m - 1);
                    if (y != 0.0) {
                        y /= x;
                        a.host(i, m - 1) = y;
                        for (int j = m; j <= n; j++) {
                            a.host(i, j) -= y * a.host(m, j);
                        }
                        for (int j = 1; j <= n; j++) {
                            a.host(j, m) += y * a.host(j, i);
                        }
                    }
                }
            }
        }
    }
}
