#pragma once
#include <cmath>
#include <matar.h>
#include "bessi1.hpp"

// Modified Bessel function K1(x) using polynomial approximation.
// Numerical Recipes (Press et al.), adapted for MATAR/Kokkos.
KOKKOS_INLINE_FUNCTION
double bessk1(double x)
{
    constexpr double p[] = {
        1.0, 0.15443144, -0.67278579, -0.18156897,
        -0.01919402, -0.00110404, -0.00004686
    };
    constexpr double q[] = {
        1.25331414, 0.23498619, -0.03655620, 0.01504268,
        -0.00780353, 0.00325614, -0.00068245
    };

    if (x <= 2.0) {
        double y = x * x / 4.0;
        return (log(x / 2.0) * bessi1(x))
             + (1.0 / x) * (p[0] + y * (p[1] + y * (p[2] + y * (p[3]
             + y * (p[4] + y * (p[5] + y * p[6]))))));
    }
    double y = 2.0 / x;
    return (exp(-x) / sqrt(x))
         * (q[0] + y * (q[1] + y * (q[2] + y * (q[3]
         + y * (q[4] + y * (q[5] + y * q[6]))))));
}
