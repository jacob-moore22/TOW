#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Safe Newton-Raphson with bisection fallback. funcd(x, f, df) sets both
// function value and derivative. Root must be bracketed in [x1,x2].
template<typename FuncD>
inline double rtsafe(FuncD funcd, double x1, double x2, double xacc)
{
    constexpr int MAXIT = 100;

    double fl, fh, df;
    funcd(x1, fl, df);
    funcd(x2, fh, df);

    if (fl * fh >= 0.0) {
        std::fprintf(stderr, "rtsafe: root must be bracketed\n");
        return 0.0;
    }

    double xl, xh;
    if (fl < 0.0) {
        xl = x1;
        xh = x2;
    } else {
        xh = x1;
        xl = x2;
        double swap = fl;
        fl = fh;
        fh = swap;
    }

    double rts   = 0.5 * (x1 + x2);
    double dxold = std::fabs(x2 - x1);
    double dx    = dxold;
    double f;
    funcd(rts, f, df);

    for (int j = 0; j < MAXIT; j++) {
        if (((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0 ||
            std::fabs(2.0 * f) > std::fabs(dxold * df)) {
            dxold = dx;
            dx    = 0.5 * (xh - xl);
            rts   = xl + dx;
            if (xl == rts) return rts;
        } else {
            dxold = dx;
            dx    = f / df;
            double temp = rts;
            rts -= dx;
            if (temp == rts) return rts;
        }
        if (std::fabs(dx) < xacc) return rts;
        funcd(rts, f, df);
        if (f < 0.0) {
            xl = rts;
            fl = f;
        } else {
            xh = rts;
            fh = f;
        }
    }

    std::fprintf(stderr, "rtsafe: exceeded maximum iterations\n");
    return rts;
}
