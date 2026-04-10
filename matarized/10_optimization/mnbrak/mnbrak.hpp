#pragma once
#include <cmath>
#include <algorithm>
#include <matar.h>

using namespace mtr;

// Bracket a minimum: given initial points ax and bx, searches downhill and
// returns a triplet (ax, bx, cx) that brackets a minimum, with corresponding
// function values (fa, fb, fc) satisfying fb <= fa and fb <= fc.
template <typename Func>
inline void mnbrak(double& ax, double& bx, double& cx,
                   double& fa, double& fb, double& fc, Func func)
{
    constexpr double GOLD   = 1.618034;
    constexpr double GLIMIT = 100.0;
    constexpr double TINY   = 1.0e-20;

    fa = func(ax);
    fb = func(bx);

    if (fb > fa) {
        std::swap(ax, bx);
        std::swap(fa, fb);
    }

    cx = bx + GOLD * (bx - ax);
    fc = func(cx);

    while (fb >= fc) {
        double r = (bx - ax) * (fb - fc);
        double q = (bx - cx) * (fb - fa);
        double u = bx - ((bx - cx) * q - (bx - ax) * r) /
                   (2.0 * std::copysign(std::max(std::fabs(q - r), TINY), q - r));
        double ulim = bx + GLIMIT * (cx - bx);
        double fu;

        if ((bx - u) * (u - cx) > 0.0) {
            fu = func(u);
            if (fu < fc) {
                ax = bx; fa = fb;
                bx = u;  fb = fu;
                return;
            } else if (fu > fb) {
                cx = u; fc = fu;
                return;
            }
            u  = cx + GOLD * (cx - bx);
            fu = func(u);
        } else if ((cx - u) * (u - ulim) > 0.0) {
            fu = func(u);
            if (fu < fc) {
                bx = cx;  cx = u;  u = cx + GOLD * (cx - bx);
                fb = fc;  fc = fu; fu = func(u);
            }
        } else if ((u - ulim) * (ulim - cx) >= 0.0) {
            u  = ulim;
            fu = func(u);
        } else {
            u  = cx + GOLD * (cx - bx);
            fu = func(u);
        }

        ax = bx;  bx = cx;  cx = u;
        fa = fb;  fb = fc;  fc = fu;
    }
}
