#pragma once
#include <cmath>
#include <matar.h>
#include "bessj1.hpp"

// Bessel function Y1(x) via rational polynomial approximation.
// Numerical Recipes (Press et al.)
// Valid for x > 0.
KOKKOS_INLINE_FUNCTION
double bessy1(double x)
{
    constexpr double p[5] = {
        1.0, 1.83105e-3, -3.516396496e-5,
        2.457520174e-6, -2.40337019e-7
    };
    constexpr double q[5] = {
        0.04687499995, -2.002690873e-4, 8.449199096e-6,
        -8.8228987e-7, 1.05787412e-7
    };
    constexpr double r[6] = {
        -4.900604943e12, 1.275274390e12, -5.153438139e10,
        7.349264551e8, -4.237922726e6, 8.511937935e3
    };
    constexpr double s[7] = {
        2.499580570e13, 4.244419664e11, 3.733650367e9,
        2.245904002e7, 1.020426050e5, 3.549632885e2, 1.0
    };

    if (x < 8.0) {
        double y = x * x;
        return x * (r[0] + y * (r[1] + y * (r[2] + y * (r[3] + y * (r[4] + y * r[5])))))
                  / (s[0] + y * (s[1] + y * (s[2] + y * (s[3] + y * (s[4] + y * (s[5] + y * s[6]))))))
             + 0.636619772 * (bessj1(x) * log(x) - 1.0 / x);
    } else {
        double z  = 8.0 / x;
        double y  = z * z;
        double xx = x - 2.356194491;
        return sqrt(0.636619772 / x)
             * (sin(xx) * (p[0] + y * (p[1] + y * (p[2] + y * (p[3] + y * p[4]))))
              + z * cos(xx) * (q[0] + y * (q[1] + y * (q[2] + y * (q[3] + y * q[4])))));
    }
}
