#pragma once
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <matar.h>

using namespace mtr;

// Brent's method. Find root of func in bracketed interval [x1,x2] to
// tolerance tol. Combines inverse quadratic interpolation, secant, and
// bisection for guaranteed convergence with superlinear speed.
template<typename Func>
inline double zbrent(Func func, double x1, double x2, double tol)
{
    constexpr int    ITMAX = 100;
    constexpr double EPS   = 3.0e-8;

    double a  = x1;
    double b  = x2;
    double fa = func(a);
    double fb = func(b);

    if (fb * fa > 0.0) {
        std::fprintf(stderr, "zbrent: root must be bracketed\n");
        return 0.0;
    }

    double c  = b;
    double fc = fb;
    double d  = b - a;
    double e  = d;

    for (int iter = 0; iter < ITMAX; iter++) {
        if (fb * fc > 0.0) {
            c  = a;
            fc = fa;
            d  = b - a;
            e  = d;
        }
        if (std::fabs(fc) < std::fabs(fb)) {
            a  = b;
            b  = c;
            c  = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        double tol1 = 2.0 * EPS * std::fabs(b) + 0.5 * tol;
        double xm   = 0.5 * (c - b);

        if (std::fabs(xm) <= tol1 || fb == 0.0) return b;

        if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
            double s = fb / fa;
            double p, q;
            if (a == c) {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                double r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0.0) q = -q;
            p = std::fabs(p);
            if (2.0 * p < std::min(3.0 * xm * q - std::fabs(tol1 * q),
                                    std::fabs(e * q))) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a  = b;
        fa = fb;
        if (std::fabs(d) > tol1) {
            b += d;
        } else {
            b += std::copysign(tol1, xm);
        }
        fb = func(b);
    }

    std::fprintf(stderr, "zbrent: exceeded maximum iterations\n");
    return b;
}
