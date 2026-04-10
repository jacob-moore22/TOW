#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "trapzd.hpp"

// Iterated trapezoidal rule to convergence.
// Returns the integral of func from a to b.
template<typename Func>
inline double qtrap(Func func, double a, double b,
                    double eps = 1.0e-6, int jmax = 20)
{
    double s, olds = -1.0e30;
    int it = 0;

    for (int j = 1; j <= jmax; j++) {
        trapzd(func, a, b, s, j, it);
        if (std::fabs(s - olds) < eps * std::fabs(olds))
            return s;
        olds = s;
    }
    std::fprintf(stderr, "qtrap: too many steps\n");
    return s;
}
