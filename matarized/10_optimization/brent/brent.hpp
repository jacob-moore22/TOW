#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Brent's method for 1-D minimization. Given a bracketing triplet (ax, bx, cx)
// with f(bx) < f(ax) and f(bx) < f(cx), finds the minimum to within tolerance
// tol. Returns the minimum value; xmin is set to the abscissa of the minimum.
template <typename Func>
inline double brent(double ax, double bx, double cx, Func f,
                    double tol, double& xmin)
{
    constexpr int    ITMAX = 100;
    constexpr double CGOLD = 0.3819660;
    constexpr double ZEPS  = 1.0e-10;

    double a = std::min(ax, cx);
    double b = std::max(ax, cx);
    double v = bx, w = bx, x = bx;
    double e = 0.0;
    double fx = f(x);
    double fv = fx, fw = fx;

    double d = 0.0;
    for (int iter = 0; iter < ITMAX; iter++) {
        double xm   = 0.5 * (a + b);
        double tol1 = tol * std::fabs(x) + ZEPS;
        double tol2 = 2.0 * tol1;

        if (std::fabs(x - xm) <= (tol2 - 0.5 * (b - a))) {
            xmin = x;
            return fx;
        }

        if (std::fabs(e) > tol1) {
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = std::fabs(q);
            double etemp = e;
            e = d;
            if (std::fabs(p) >= std::fabs(0.5 * q * etemp) ||
                p <= q * (a - x) || p >= q * (b - x)) {
                e = (x >= xm) ? a - x : b - x;
                d = CGOLD * e;
            } else {
                d = p / q;
                double u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = std::copysign(tol1, xm - x);
            }
        } else {
            e = (x >= xm) ? a - x : b - x;
            d = CGOLD * e;
        }

        double u = (std::fabs(d) >= tol1) ? x + d : x + std::copysign(tol1, d);
        double fu = f(u);

        if (fu <= fx) {
            if (u >= x) a = x; else b = x;
            v = w;  fv = fw;
            w = x;  fw = fx;
            x = u;  fx = fu;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w;  fv = fw;
                w = u;  fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u;  fv = fu;
            }
        }
    }
    std::fprintf(stderr, "brent: exceeded maximum iterations\n");
    xmin = x;
    return fx;
}
