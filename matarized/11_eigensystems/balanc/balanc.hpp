// Balance a nonsymmetric matrix by row and column scaling to reduce its norm
// (Numerical Recipes BALANC).
// Replaces a[1..n][1..n] by a balanced matrix with identical eigenvalues.
// RADIX should be the machine's floating-point radix.

#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

inline void balanc(DFMatrixKokkos<double>& a, int n)
{
    constexpr double RADIX = 2.0;
    constexpr double SQRDX = RADIX * RADIX;

    int last = 0;
    while (last == 0) {
        last = 1;
        for (int i = 1; i <= n; i++) {
            double c = 0.0;
            double r = 0.0;
            for (int j = 1; j <= n; j++) {
                if (j != i) {
                    c += std::fabs(a.host(j, i));
                    r += std::fabs(a.host(i, j));
                }
            }
            if (c != 0.0 && r != 0.0) {
                double g = r / RADIX;
                double f = 1.0;
                double s = c + r;
                while (c < g) {
                    f *= RADIX;
                    c *= SQRDX;
                }
                g = r * RADIX;
                while (c > g) {
                    f /= RADIX;
                    c /= SQRDX;
                }
                if ((c + r) / f < 0.95 * s) {
                    last = 0;
                    g = 1.0 / f;
                    for (int j = 1; j <= n; j++) {
                        a.host(i, j) *= g;
                    }
                    for (int j = 1; j <= n; j++) {
                        a.host(j, i) *= f;
                    }
                }
            }
        }
    }
}
