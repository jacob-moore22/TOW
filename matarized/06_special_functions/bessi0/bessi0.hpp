#pragma once
#include <cmath>
#include <matar.h>

// Modified Bessel function I0(x) using polynomial approximation.
// Numerical Recipes (Press et al.), adapted for MATAR/Kokkos.
KOKKOS_INLINE_FUNCTION
double bessi0(double x)
{
    constexpr double p[] = {
        1.0, 3.5156229, 3.0899424, 1.2067492,
        0.2659732, 0.0360768, 0.0045813
    };
    constexpr double q[] = {
        0.39894228, 0.01328592, 0.00225319, -0.00157565,
        0.00916281, -0.02057706, 0.02635537, -0.01647633,
        0.00392377
    };

    double ax = fabs(x);
    if (ax < 3.75) {
        double y = (x / 3.75) * (x / 3.75);
        return p[0] + y * (p[1] + y * (p[2] + y * (p[3]
             + y * (p[4] + y * (p[5] + y * p[6])))));
    }
    double y = 3.75 / ax;
    return (exp(ax) / sqrt(ax))
         * (q[0] + y * (q[1] + y * (q[2] + y * (q[3]
         + y * (q[4] + y * (q[5] + y * (q[6] + y * (q[7]
         + y * q[8]))))))));
}
