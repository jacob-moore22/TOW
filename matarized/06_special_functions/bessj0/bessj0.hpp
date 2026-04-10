#pragma once
#include <cmath>
#include <matar.h>

// Bessel function J0(x) via rational polynomial approximation.
// Numerical Recipes (Press et al.)
KOKKOS_INLINE_FUNCTION
double bessj0(double x)
{
    constexpr double p[5] = {
        1.0, -1.098628627e-3, 2.734510407e-5,
        -2.073370639e-6, 2.093887211e-7
    };
    constexpr double q[5] = {
        -1.562499995e-2, 1.430488765e-4, -6.911147651e-6,
        7.621095161e-7, -9.34945152e-8
    };
    constexpr double r[6] = {
        57568490574.0, -13362590354.0, 651619640.7,
        -11214424.18, 77392.33017, -184.9052456
    };
    constexpr double s[6] = {
        57568490411.0, 1029532985.0, 9494680.718,
        59272.64853, 267.8532712, 1.0
    };

    if (fabs(x) < 8.0) {
        double y = x * x;
        return (r[0] + y * (r[1] + y * (r[2] + y * (r[3] + y * (r[4] + y * r[5])))))
             / (s[0] + y * (s[1] + y * (s[2] + y * (s[3] + y * (s[4] + y * s[5])))));
    } else {
        double ax = fabs(x);
        double z  = 8.0 / ax;
        double y  = z * z;
        double xx = ax - 0.785398164;
        return sqrt(0.636619772 / ax)
             * (cos(xx) * (p[0] + y * (p[1] + y * (p[2] + y * (p[3] + y * p[4]))))
              - z * sin(xx) * (q[0] + y * (q[1] + y * (q[2] + y * (q[3] + y * q[4])))));
    }
}
