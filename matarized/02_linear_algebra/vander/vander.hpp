// Vandermonde system solver (Numerical Recipes VANDER).
//
// Solves the Vandermonde linear system:
//   sum_{i=1}^{n} x(k)^(i-1) * w(i) = q(k)  for k=1..n
//
// x(1..n) -- the abscissas (distinct points)
// w(1..n) -- solution (output): the coefficients
// q(1..n) -- right-hand side values
// n       -- system size

#pragma once
#include <cstdio>
#include <matar.h>

using namespace mtr;

inline void vander(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& w,
                   DFMatrixKokkos<double>& q, int n)
{
    if (n == 1) {
        w.host(1) = q.host(1);
        return;
    }

    double c[n + 1];
    for (int i = 1; i <= n; i++) c[i] = 0.0;

    c[n] = -x.host(1);
    for (int i = 2; i <= n; i++) {
        double xx = -x.host(i);
        for (int j = n + 1 - i; j <= n - 1; j++)
            c[j] += xx * c[j + 1];
        c[n] += xx;
    }

    for (int i = 1; i <= n; i++) {
        double xx = x.host(i);
        double t = 1.0;
        double b = 1.0;
        double s = q.host(n);
        int k = n;
        for (int j = 2; j <= n; j++) {
            int k1 = k - 1;
            b = c[k] + xx * b;
            s += q.host(k1) * b;
            t = xx * t + b;
            k = k1;
        }
        w.host(i) = s / t;
    }
}
