#pragma once
#include <cmath>
#include <matar.h>

// Bessel function J1(x) via rational polynomial approximation.
// Numerical Recipes (Press et al.)
KOKKOS_INLINE_FUNCTION
double bessj1(double x)
{
    constexpr double r[6] = {
        72362614232.0, -7895059235.0, 242396853.1,
        -2972611.439, 15704.48260, -30.16036606
    };
    constexpr double s[6] = {
        144725228442.0, 2300535178.0, 18583304.74,
        99447.43394, 376.9991397, 1.0
    };
    constexpr double p[5] = {
        1.0, 1.83105e-3, -3.516396496e-5,
        2.457520174e-6, -2.40337019e-7
    };
    constexpr double q[5] = {
        0.04687499995, -2.002690873e-4, 8.449199096e-6,
        -8.8228987e-7, 1.05787412e-7
    };

    if (fabs(x) < 8.0) {
        double y = x * x;
        return x * (r[0] + y * (r[1] + y * (r[2] + y * (r[3] + y * (r[4] + y * r[5])))))
                  / (s[0] + y * (s[1] + y * (s[2] + y * (s[3] + y * (s[4] + y * s[5])))));
    } else {
        double ax = fabs(x);
        double z  = 8.0 / ax;
        double y  = z * z;
        double xx = ax - 2.356194491;
        double sign = (x >= 0.0) ? 1.0 : -1.0;
        return sqrt(0.636619772 / ax)
             * (cos(xx) * (p[0] + y * (p[1] + y * (p[2] + y * (p[3] + y * p[4]))))
              - z * sin(xx) * (q[0] + y * (q[1] + y * (q[2] + y * (q[3] + y * q[4])))))
             * sign;
    }
}
