// Polynomial coefficients from data (Numerical Recipes POLCOE).
// Given arrays x[1..n] and y[1..n], returns polynomial coefficients
// cof[1..n] such that y = cof(1) + cof(2)*x + ... + cof(n)*x^(n-1).

#pragma once
#include <matar.h>

using namespace mtr;

inline void polcoe(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                   int n, DFMatrixKokkos<double>& cof)
{
    constexpr int NMAX = 15;
    double s[NMAX + 1];

    for (int i = 1; i <= n; i++) {
        s[i]          = 0.0;
        cof.host(i)   = 0.0;
    }
    s[n] = -x.host(1);
    for (int i = 2; i <= n; i++) {
        for (int j = n + 1 - i; j <= n - 1; j++) {
            s[j] = s[j] - x.host(i) * s[j + 1];
        }
        s[n] = s[n] - x.host(i);
    }
    for (int j = 1; j <= n; j++) {
        double phi = static_cast<double>(n);
        for (int k = n - 1; k >= 1; k--) {
            phi = k * s[k + 1] + x.host(j) * phi;
        }
        double ff = y.host(j) / phi;
        double b  = 1.0;
        for (int k = n; k >= 1; k--) {
            cof.host(k) = cof.host(k) + b * ff;
            b = s[k] + x.host(j) * b;
        }
    }
}
