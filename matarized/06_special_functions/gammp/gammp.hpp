#pragma once
#include <cmath>
#include <matar.h>
#include "gser.hpp"
#include "gcf.hpp"

// Incomplete gamma function P(a,x) = gamma(a,x) / Gamma(a).
KOKKOS_INLINE_FUNCTION
double gammp(double a, double x)
{
    double gln;
    if (x < a + 1.0) {
        double gamser;
        gser(gamser, a, x, gln);
        return gamser;
    } else {
        double gammcf;
        gcf(gammcf, a, x, gln);
        return 1.0 - gammcf;
    }
}
