#pragma once
#include <matar.h>

using namespace mtr;

// Shift polynomial coefficients d(1..n) from the Chebyshev variable
// y in [-1,1] to the original variable x in [a,b].
// After this call, d(1) + d(2)*x + d(3)*x^2 + ... approximates f(x).
inline void pcshft(double a, double b, DFMatrixKokkos<double>& d, int n)
{
    double cnst = 2.0 / (b - a);
    double fac  = cnst;

    for (int j = 2; j <= n; j++) {
        d(j) *= fac;
        fac  *= cnst;
    }

    cnst = 0.5 * (a + b);
    for (int j = 1; j <= n - 1; j++) {
        for (int k = n - 1; k >= j; k--)
            d(k) -= cnst * d(k + 1);
    }
}
