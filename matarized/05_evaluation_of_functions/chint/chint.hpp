#pragma once
#include <matar.h>

using namespace mtr;

// Integral of Chebyshev series: given c(1..n) on [a,b], produce cint(1..n).
// Evaluating the resulting series gives the integral from a to x of the original.
inline void chint(double a, double b,
                  const DFMatrixKokkos<double>& c,
                  DFMatrixKokkos<double>& cint, int n)
{
    double con = 0.25 * (b - a);
    double sum = 0.0;
    double fac = 1.0;

    for (int j = 2; j <= n - 1; j++) {
        cint(j) = con * (c(j - 1) - c(j + 1)) / (j - 1);
        sum += fac * cint(j);
        fac = -fac;
    }

    cint(n) = con * c(n - 1) / (n - 1);
    sum += fac * cint(n);
    cint(1) = 2.0 * sum;
}
