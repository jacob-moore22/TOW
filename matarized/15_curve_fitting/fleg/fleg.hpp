// Legendre polynomial basis functions (Numerical Recipes FLEG).
//
// Evaluates pl(1..nl) = {P_0(x), P_1(x), ..., P_{nl-1}(x)} using
// the standard recurrence relation. For use with lfit/svdfit.

#pragma once
#include <matar.h>

using namespace mtr;

inline void fleg(double x, DFMatrixKokkos<double>& pl, int nl)
{
    pl.host(1) = 1.0;
    pl.host(2) = x;
    if (nl > 2) {
        double twox = 2.0 * x;
        double f2 = x;
        double d = 1.0;
        for (int j = 3; j <= nl; j++) {
            double f1 = d;
            f2 += twox;
            d += 1.0;
            pl.host(j) = (f2 * pl.host(j - 1) - f1 * pl.host(j - 2)) / d;
        }
    }
}
