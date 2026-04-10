#pragma once
#include <matar.h>

// 10-point Gauss-Legendre quadrature.
// Returns the integral of func from a to b.
template<typename Func>
KOKKOS_INLINE_FUNCTION
double qgaus(Func func, double a, double b)
{
    constexpr double x[5] = {
        0.1488743389, 0.4333953941, 0.6794095682,
        0.8650633666, 0.9739065285
    };
    constexpr double w[5] = {
        0.2955242247, 0.2692667193, 0.2190863625,
        0.1494513491, 0.0666713443
    };

    double xm = 0.5 * (b + a);
    double xr = 0.5 * (b - a);
    double ss = 0.0;
    for (int j = 0; j < 5; j++) {
        double dx = xr * x[j];
        ss += w[j] * (func(xm + dx) + func(xm - dx));
    }
    return xr * ss;
}
