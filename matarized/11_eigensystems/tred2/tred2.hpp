// Householder reduction of a real symmetric matrix to tridiagonal form
// (Numerical Recipes TRED2).
// On input, a[1..n][1..n] is the symmetric matrix. On output, a is replaced
// by the orthogonal transformation matrix Q, d[1..n] contains the diagonal
// elements of the tridiagonal matrix, and e[1..n] the off-diagonal elements
// (with e[1]=0).

#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

inline void tred2(DFMatrixKokkos<double>& a, int n,
                  DFMatrixKokkos<double>& d, DFMatrixKokkos<double>& e)
{
    if (n > 1) {
        for (int i = n; i >= 2; i--) {
            int l = i - 1;
            double h = 0.0;
            double scale = 0.0;
            if (l > 1) {
                for (int k = 1; k <= l; k++) {
                    scale += std::fabs(a.host(i, k));
                }
                if (scale == 0.0) {
                    e.host(i) = a.host(i, l);
                } else {
                    for (int k = 1; k <= l; k++) {
                        a.host(i, k) /= scale;
                        h += a.host(i, k) * a.host(i, k);
                    }
                    double f = a.host(i, l);
                    double g = -std::copysign(std::sqrt(h), f);
                    e.host(i) = scale * g;
                    h -= f * g;
                    a.host(i, l) = f - g;
                    f = 0.0;
                    for (int j = 1; j <= l; j++) {
                        a.host(j, i) = a.host(i, j) / h;
                        g = 0.0;
                        for (int k = 1; k <= j; k++) {
                            g += a.host(j, k) * a.host(i, k);
                        }
                        if (l > j) {
                            for (int k = j + 1; k <= l; k++) {
                                g += a.host(k, j) * a.host(i, k);
                            }
                        }
                        e.host(j) = g / h;
                        f += e.host(j) * a.host(i, j);
                    }
                    double hh = f / (h + h);
                    for (int j = 1; j <= l; j++) {
                        f = a.host(i, j);
                        g = e.host(j) - hh * f;
                        e.host(j) = g;
                        for (int k = 1; k <= j; k++) {
                            a.host(j, k) -= f * e.host(k) + g * a.host(i, k);
                        }
                    }
                }
            } else {
                e.host(i) = a.host(i, l);
            }
            d.host(i) = h;
        }
    }
    d.host(1) = 0.0;
    e.host(1) = 0.0;
    for (int i = 1; i <= n; i++) {
        int l = i - 1;
        if (d.host(i) != 0.0) {
            for (int j = 1; j <= l; j++) {
                double g = 0.0;
                for (int k = 1; k <= l; k++) {
                    g += a.host(i, k) * a.host(k, j);
                }
                for (int k = 1; k <= l; k++) {
                    a.host(k, j) -= g * a.host(k, i);
                }
            }
        }
        d.host(i) = a.host(i, i);
        a.host(i, i) = 1.0;
        if (l >= 1) {
            for (int j = 1; j <= l; j++) {
                a.host(i, j) = 0.0;
                a.host(j, i) = 0.0;
            }
        }
    }
}
