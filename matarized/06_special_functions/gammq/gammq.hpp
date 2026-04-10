#pragma once
#include <cmath>
#include <matar.h>
#include "gser.hpp"
#include "gcf.hpp"

// Complementary incomplete gamma function Q(a,x) = 1 - P(a,x).
KOKKOS_INLINE_FUNCTION
double gammq(double a, double x)
{
    double gln;
    if (x < a + 1.0) {
        double gamser;
        gser(gamser, a, x, gln);
        return 1.0 - gamser;
    } else {
        double gammcf;
        gcf(gammcf, a, x, gln);
        return gammcf;
    }
}
