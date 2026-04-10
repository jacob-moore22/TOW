#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Bracket a root of func. Given initial range [x1,x2], expands the interval
// outward until func(x1)*func(x2) < 0 (i.e. a sign change is found).
// Returns true on success, false if no bracket found after NTRY attempts.
template<typename Func>
inline bool zbrac(Func func, double& x1, double& x2)
{
    constexpr double FACTOR = 1.6;
    constexpr int    NTRY   = 50;

    if (x1 == x2) {
        std::fprintf(stderr, "zbrac: initial range has zero width\n");
        return false;
    }

    double f1 = func(x1);
    double f2 = func(x2);

    for (int j = 0; j < NTRY; j++) {
        if (f1 * f2 < 0.0) return true;
        if (std::fabs(f1) < std::fabs(f2)) {
            x1 += FACTOR * (x1 - x2);
            f1 = func(x1);
        } else {
            x2 += FACTOR * (x2 - x1);
            f2 = func(x2);
        }
    }
    return false;
}
