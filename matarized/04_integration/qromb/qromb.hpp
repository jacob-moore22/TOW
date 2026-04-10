#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "trapzd.hpp"
#include "polint.hpp"

// Romberg integration using the extended trapezoidal rule and
// polynomial extrapolation (polint).  Returns the integral of
// func from a to b.
template<typename Func>
inline double qromb(Func func, double a, double b,
                    double eps = 1.0e-6, int jmax = 20)
{
    constexpr int K  = 5;
    constexpr int KM = K - 1;
    constexpr int JMAXP = 21;

    double s[JMAXP], h[JMAXP];
    h[0] = 1.0;
    int it = 0;
    double ss = 0.0;

    for (int j = 0; j < jmax; j++) {
        trapzd(func, a, b, s[j], j + 1, it);
        if (j + 1 >= K) {
            double dss;
            polint(&h[j + 1 - K], &s[j + 1 - K], K, 0.0, ss, dss);
            if (std::fabs(dss) < eps * std::fabs(ss))
                return ss;
        }
        s[j + 1] = s[j];
        h[j + 1] = 0.25 * h[j];
    }
    std::fprintf(stderr, "qromb: too many steps\n");
    return ss;
}
