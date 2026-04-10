// QL algorithm with implicit shifts for eigenvalues of a real symmetric
// tridiagonal matrix (Numerical Recipes TQLI).
// On input, d[1..n] contains the diagonal and e[1..n] the off-diagonal
// (e[1] arbitrary) of the tridiagonal matrix. z[1..n][1..n] should be
// the identity (for eigenvalues of the tridiagonal matrix) or the output
// of tred2 (for eigenvalues of the original matrix). On output, d contains
// eigenvalues and columns of z are the eigenvectors.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

inline void tqli(DFMatrixKokkos<double>& d, DFMatrixKokkos<double>& e,
                 int n, DFMatrixKokkos<double>& z)
{
    if (n > 1) {
        for (int i = 2; i <= n; i++) {
            e.host(i - 1) = e.host(i);
        }
        e.host(n) = 0.0;

        for (int l = 1; l <= n; l++) {
            int iter = 0;
            int m;
            while (true) {
                for (m = l; m <= n - 1; m++) {
                    double dd = std::fabs(d.host(m)) + std::fabs(d.host(m + 1));
                    if (std::fabs(e.host(m)) + dd == dd) break;
                }
                if (m == l) break;

                if (iter == 30) {
                    std::fprintf(stderr, "tqli: too many iterations\n");
                    return;
                }
                iter++;

                double g = (d.host(l + 1) - d.host(l)) / (2.0 * e.host(l));
                double r = std::sqrt(g * g + 1.0);
                g = d.host(m) - d.host(l) + e.host(l) / (g + std::copysign(r, g));
                double s = 1.0;
                double c = 1.0;
                double p = 0.0;

                for (int i = m - 1; i >= l; i--) {
                    double f = s * e.host(i);
                    double b = c * e.host(i);
                    if (std::fabs(f) >= std::fabs(g)) {
                        c = g / f;
                        r = std::sqrt(c * c + 1.0);
                        e.host(i + 1) = f * r;
                        s = 1.0 / r;
                        c *= s;
                    } else {
                        s = f / g;
                        r = std::sqrt(s * s + 1.0);
                        e.host(i + 1) = g * r;
                        c = 1.0 / r;
                        s *= c;
                    }
                    g = d.host(i + 1) - p;
                    r = (d.host(i) - g) * s + 2.0 * c * b;
                    p = s * r;
                    d.host(i + 1) = g + p;
                    g = c * r - b;

                    for (int k = 1; k <= n; k++) {
                        f = z.host(k, i + 1);
                        z.host(k, i + 1) = s * z.host(k, i) + c * f;
                        z.host(k, i) = c * z.host(k, i) - s * f;
                    }
                }
                d.host(l) -= p;
                e.host(l) = g;
                e.host(m) = 0.0;
            }
        }
    }
}
