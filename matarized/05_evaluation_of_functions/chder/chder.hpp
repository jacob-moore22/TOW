#pragma once
#include <matar.h>

using namespace mtr;

// Derivative of Chebyshev series: given c(1..n) on [a,b], produce cder(1..n).
inline void chder(double a, double b,
                  const DFMatrixKokkos<double>& c,
                  DFMatrixKokkos<double>& cder, int n)
{
    cder(n)     = 0.0;
    cder(n - 1) = 2.0 * (n - 1) * c(n);

    if (n >= 3) {
        for (int j = n - 2; j >= 1; j--)
            cder(j) = cder(j + 2) + 2.0 * j * c(j + 1);
    }

    double con = 2.0 / (b - a);
    for (int j = 1; j <= n; j++)
        cder(j) *= con;
}
