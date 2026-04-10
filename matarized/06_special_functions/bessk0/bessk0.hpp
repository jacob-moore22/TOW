#pragma once
#include <cmath>
#include <matar.h>
#include "bessi0.hpp"

// Modified Bessel function K0(x) using polynomial approximation.
// Numerical Recipes (Press et al.), adapted for MATAR/Kokkos.
KOKKOS_INLINE_FUNCTION
double bessk0(double x)
{
    constexpr double p[] = {
        -0.57721566, 0.42278420, 0.23069756, 0.03488590,
         0.00262698, 0.00010750, 0.0000074
    };
    constexpr double q[] = {
        1.25331414, -0.07832358, 0.02189568, -0.01062446,
        0.00587872, -0.00251540, 0.00053208
    };

    if (x <= 2.0) {
        double y = x * x / 4.0;
        return (-log(x / 2.0) * bessi0(x))
             + (p[0] + y * (p[1] + y * (p[2] + y * (p[3]
             + y * (p[4] + y * (p[5] + y * p[6]))))));
    }
    double y = 2.0 / x;
    return (exp(-x) / sqrt(x))
         * (q[0] + y * (q[1] + y * (q[2] + y * (q[3]
         + y * (q[4] + y * (q[5] + y * q[6]))))));
}
