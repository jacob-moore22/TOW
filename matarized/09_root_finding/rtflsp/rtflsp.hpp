#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// False position (regula falsi) method. Find root of func in bracketed
// interval [x1,x2] to accuracy xacc.
template<typename Func>
inline double rtflsp(Func func, double x1, double x2, double xacc)
{
    constexpr int MAXIT = 30;

    double fl = func(x1);
    double fh = func(x2);

    if (fl * fh > 0.0) {
        std::fprintf(stderr, "rtflsp: root must be bracketed\n");
        return 0.0;
    }

    double xl, xh;
    if (fl < 0.0) {
        xl = x1;
        xh = x2;
    } else {
        xl = x2;
        xh = x1;
        double swap = fl;
        fl = fh;
        fh = swap;
    }

    double dx = xh - xl;

    for (int j = 0; j < MAXIT; j++) {
        double rt = xl + dx * fl / (fl - fh);
        double f  = func(rt);
        double del;
        if (f < 0.0) {
            del = xl - rt;
            xl  = rt;
            fl  = f;
        } else {
            del = xh - rt;
            xh  = rt;
            fh  = f;
        }
        dx = xh - xl;
        if (std::fabs(del) < xacc || f == 0.0) return rt;
    }

    std::fprintf(stderr, "rtflsp: exceeded maximum iterations\n");
    return xl + dx * fl / (fl - fh);
}
