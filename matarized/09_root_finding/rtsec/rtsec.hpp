#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Secant method. Find root of func given two initial guesses x1, x2
// (not necessarily bracketing) to accuracy xacc.
template<typename Func>
inline double rtsec(Func func, double x1, double x2, double xacc)
{
    constexpr int MAXIT = 30;

    double fl = func(x1);
    double f  = func(x2);

    double rts, xl;
    if (std::fabs(fl) < std::fabs(f)) {
        rts = x1;
        xl  = x2;
        double swap = fl;
        fl = f;
        f  = swap;
    } else {
        xl  = x1;
        rts = x2;
    }

    for (int j = 0; j < MAXIT; j++) {
        double dx = (xl - rts) * f / (f - fl);
        xl  = rts;
        fl  = f;
        rts += dx;
        f   = func(rts);
        if (std::fabs(dx) < xacc || f == 0.0) return rts;
    }

    std::fprintf(stderr, "rtsec: exceeded maximum iterations\n");
    return rts;
}
