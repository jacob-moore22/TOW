#pragma once
#include <cmath>
#include <matar.h>

// Modified Bessel function I1(x) using polynomial approximation.
// Numerical Recipes (Press et al.), adapted for MATAR/Kokkos.
KOKKOS_INLINE_FUNCTION
double bessi1(double x)
{
    constexpr double p[] = {
        0.5, 0.87890594, 0.51498869, 0.15084934,
        0.02658733, 0.00301532, 0.00032411
    };
    constexpr double q[] = {
        0.39894228, -0.03988024, -0.00362018, 0.00163801,
        -0.01031555, 0.02282967, -0.02895312, 0.01787654,
        -0.00420059
    };

    double ax = fabs(x);
    if (ax < 3.75) {
        double y = (x / 3.75) * (x / 3.75);
        return x * (p[0] + y * (p[1] + y * (p[2] + y * (p[3]
             + y * (p[4] + y * (p[5] + y * p[6]))))));
    }
    double y = 3.75 / ax;
    double ans = (exp(ax) / sqrt(ax))
              * (q[0] + y * (q[1] + y * (q[2] + y * (q[3]
              + y * (q[4] + y * (q[5] + y * (q[6] + y * (q[7]
              + y * q[8]))))))));
    return (x < 0.0) ? -ans : ans;
}
