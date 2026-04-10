#pragma once
#include <cmath>
#include <matar.h>
#include "gammln.hpp"

// Incomplete gamma function Q(a,x) evaluated by its continued-fraction representation.
// Returns gammcf = Q(a,x) and gln = ln(Gamma(a)).
KOKKOS_INLINE_FUNCTION
void gcf(double& gammcf, double a, double x, double& gln)
{
    constexpr int ITMAX = 100;
    constexpr double EPS = 3.0e-7;

    gln = gammln(a);

    double gold = 0.0;
    double a0   = 1.0;
    double a1   = x;
    double b0   = 0.0;
    double b1   = 1.0;
    double fac  = 1.0;
    double g    = 0.0;

    for (int n = 1; n <= ITMAX; n++) {
        double an  = double(n);
        double ana = an - a;
        a0 = (a1 + a0 * ana) * fac;
        b0 = (b1 + b0 * ana) * fac;
        double anf = an * fac;
        a1 = x * a0 + anf * a1;
        b1 = x * b0 + anf * b1;
        if (a1 != 0.0) {
            fac = 1.0 / a1;
            g = b1 * fac;
            if (fabs((g - gold) / g) < EPS) {
                gammcf = exp(-x + a * log(x) - gln) * g;
                return;
            }
            gold = g;
        }
    }
    gammcf = exp(-x + a * log(x) - gln) * g;
}
