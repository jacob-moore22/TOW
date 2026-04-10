#pragma once
#include <cmath>
#include <matar.h>
#include "gammln.hpp"

// Incomplete gamma function P(a,x) evaluated by its series representation.
// Returns gamser = P(a,x) and gln = ln(Gamma(a)).
KOKKOS_INLINE_FUNCTION
void gser(double& gamser, double a, double x, double& gln)
{
    constexpr int ITMAX = 100;
    constexpr double EPS = 3.0e-7;

    gln = gammln(a);

    if (x <= 0.0) {
        gamser = 0.0;
        return;
    }

    double ap  = a;
    double sum = 1.0 / a;
    double del = sum;

    for (int n = 1; n <= ITMAX; n++) {
        ap  += 1.0;
        del *= x / ap;
        sum += del;
        if (fabs(del) < fabs(sum) * EPS) {
            gamser = sum * exp(-x + a * log(x) - gln);
            return;
        }
    }
    gamser = sum * exp(-x + a * log(x) - gln);
}
