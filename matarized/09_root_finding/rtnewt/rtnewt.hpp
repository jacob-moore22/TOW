#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Newton-Raphson root finding. funcd(x, f, df) sets the function value f
// and derivative df at x. Searches within [x1,x2] to accuracy xacc.
template<typename FuncD>
inline double rtnewt(FuncD funcd, double x1, double x2, double xacc)
{
    constexpr int JMAX = 20;

    double rtn = 0.5 * (x1 + x2);

    for (int j = 0; j < JMAX; j++) {
        double f, df;
        funcd(rtn, f, df);
        double dx = f / df;
        rtn -= dx;
        if ((x1 - rtn) * (rtn - x2) < 0.0) {
            std::fprintf(stderr, "rtnewt: jumped out of brackets\n");
            return rtn;
        }
        if (std::fabs(dx) < xacc) return rtn;
    }

    std::fprintf(stderr, "rtnewt: exceeded maximum iterations\n");
    return rtn;
}
