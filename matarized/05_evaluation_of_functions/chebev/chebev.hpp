#pragma once
#include <matar.h>

using namespace mtr;

// Evaluate Chebyshev series c(1..m) at point x on interval [a,b].
KOKKOS_INLINE_FUNCTION
double chebev(double a, double b, const DFMatrixKokkos<double>& c, int m, double x)
{
    double d  = 0.0;
    double dd = 0.0;
    double y  = (2.0 * x - a - b) / (b - a);
    double y2 = 2.0 * y;

    for (int j = m; j >= 2; j--) {
        double sv = d;
        d  = y2 * d - dd + c(j);
        dd = sv;
    }
    return y * d - dd + 0.5 * c(1);
}
