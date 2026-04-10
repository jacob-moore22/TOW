#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Bisection method. Find root of func in bracketed interval [x1,x2]
// to accuracy xacc. func(x1) and func(x2) must have opposite signs.
template<typename Func>
inline double rtbis(Func func, double x1, double x2, double xacc)
{
    constexpr int JMAX = 40;

    double fmid = func(x2);
    double f    = func(x1);

    if (f * fmid >= 0.0) {
        std::fprintf(stderr, "rtbis: root must be bracketed\n");
        return 0.0;
    }

    double rtb, dx;
    if (f < 0.0) {
        rtb = x1;
        dx  = x2 - x1;
    } else {
        rtb = x2;
        dx  = x1 - x2;
    }

    for (int j = 0; j < JMAX; j++) {
        dx  *= 0.5;
        double xmid = rtb + dx;
        fmid = func(xmid);
        if (fmid <= 0.0) rtb = xmid;
        if (std::fabs(dx) < xacc || fmid == 0.0) return rtb;
    }

    std::fprintf(stderr, "rtbis: too many bisections\n");
    return rtb;
}
