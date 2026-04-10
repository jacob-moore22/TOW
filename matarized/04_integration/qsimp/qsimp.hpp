#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "trapzd.hpp"

// Simpson's rule to convergence using the trapezoidal refinement.
// Returns the integral of func from a to b.
template<typename Func>
inline double qsimp(Func func, double a, double b,
                    double eps = 1.0e-6, int jmax = 20)
{
    double s, st, ost = -1.0e30, os = -1.0e30;
    int it = 0;

    for (int j = 1; j <= jmax; j++) {
        trapzd(func, a, b, st, j, it);
        s = (4.0 * st - ost) / 3.0;
        if (std::fabs(s - os) < eps * std::fabs(os))
            return s;
        os  = s;
        ost = st;
    }
    std::fprintf(stderr, "qsimp: too many steps\n");
    return s;
}
