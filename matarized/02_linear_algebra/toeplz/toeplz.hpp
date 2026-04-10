// Toeplitz system solver using Levinson recursion (Numerical Recipes TOEPLZ).
//
// Solves the Toeplitz system:  sum_{j=1}^{n} R(n+i-j) * x(j) = y(i)
//
// r(1..2n-1) -- row of the Toeplitz matrix: r(n) is the diagonal,
//               r(n-1)..r(1) are sub-diagonals, r(n+1)..r(2n-1) super-diags.
// x(1..n)    -- solution (output)
// y(1..n)    -- right-hand side
// n          -- system size

#pragma once
#include <cstdio>
#include <matar.h>

using namespace mtr;

inline void toeplz(DFMatrixKokkos<double>& r, DFMatrixKokkos<double>& x,
                   DFMatrixKokkos<double>& y, int n)
{
    double g[n], h[n];

    if (r.host(n) == 0.0) {
        std::fprintf(stderr, "toeplz: singular principal minor\n");
        return;
    }

    x.host(1) = y.host(1) / r.host(n);
    if (n == 1) return;

    g[0] = r.host(n - 1) / r.host(n);
    h[0] = r.host(n + 1) / r.host(n);

    for (int m = 1; m < n; m++) {
        int m1 = m + 1;

        double sxn = -y.host(m1);
        double sd  = -r.host(n);
        for (int j = 1; j <= m; j++) {
            sxn += r.host(n + m1 - j) * x.host(j);
            sd  += r.host(n + m1 - j) * g[m - j];
        }
        if (sd == 0.0) {
            std::fprintf(stderr, "toeplz: singular principal minor\n");
            return;
        }
        x.host(m1) = sxn / sd;
        for (int j = 1; j <= m; j++)
            x.host(j) -= x.host(m1) * g[m - j];

        if (m1 == n) return;

        double sgn = -r.host(n - m1);
        double shn = -r.host(n + m1);
        double sgd = -r.host(n);
        for (int j = 1; j <= m; j++) {
            sgn += r.host(n + j - m1) * g[j - 1];
            shn += r.host(n + m1 - j) * h[j - 1];
            sgd += r.host(n + j - m1) * h[m - j];
        }
        if (sgd == 0.0) {
            std::fprintf(stderr, "toeplz: singular principal minor\n");
            return;
        }
        g[m] = sgn / sgd;
        h[m] = shn / sd;

        int k = m;
        int m2 = (m + 1) / 2;
        double pp = g[m];
        double qq = h[m];
        for (int j = 0; j < m2; j++) {
            double pt1 = g[j];
            double pt2 = g[k - 1];
            double qt1 = h[j];
            double qt2 = h[k - 1];
            g[j]     = pt1 - pp * qt2;
            g[k - 1] = pt2 - pp * qt1;
            h[j]     = qt1 - qq * pt2;
            h[k - 1] = qt2 - qq * pt1;
            k--;
        }
    }

    std::fprintf(stderr, "toeplz: should not reach here\n");
}
