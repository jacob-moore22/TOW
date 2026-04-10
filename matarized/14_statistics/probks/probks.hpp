#pragma once
#include <cmath>
#include <matar.h>

// Kolmogorov-Smirnov probability function Q_KS(lambda).
KOKKOS_INLINE_FUNCTION
double probks(double alam)
{
    constexpr double EPS1 = 0.001;
    constexpr double EPS2 = 1.0e-8;

    double a2  = -2.0 * alam * alam;
    double fac = 2.0;
    double prob = 0.0;
    double termbf = 0.0;

    for (int j = 1; j <= 100; j++) {
        double term = fac * exp(a2 * j * j);
        prob += term;
        if (fabs(term) < EPS1 * termbf || fabs(term) < EPS2 * prob)
            return prob;
        fac = -fac;
        termbf = fabs(term);
    }
    return 1.0;
}
