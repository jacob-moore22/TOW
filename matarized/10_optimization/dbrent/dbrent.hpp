#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Brent's method with derivative information. Given a bracketing triplet
// (ax, bx, cx), function f, and derivative df, finds the minimum to within
// tolerance tol. Returns the minimum value; xmin receives the location.
template <typename Func, typename DFunc>
inline double dbrent(double ax, double bx, double cx, Func f, DFunc df,
                     double tol, double& xmin)
{
    constexpr int    ITMAX = 100;
    constexpr double ZEPS  = 1.0e-10;

    double a = std::min(ax, cx);
    double b = std::max(ax, cx);
    double v = bx, w = bx, x = bx;
    double e = 0.0;
    double fx = f(x), fv = fx, fw = fx;
    double dx = df(x), dv = dx, dw = dx;

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
            double d1 = 2.0 * (b - a);
            double d2 = d1;
            if (dw != dx) d1 = (w - x) * dx / (dx - dw);
            if (dv != dx) d2 = (v - x) * dx / (dx - dv);
            double u1 = x + d1;
            double u2 = x + d2;
            bool ok1 = ((a - u1) * (u1 - b) > 0.0) && (dx * d1 <= 0.0);
            bool ok2 = ((a - u2) * (u2 - b) > 0.0) && (dx * d2 <= 0.0);
            double olde = e;
            e = d;

            if (!ok1 && !ok2) {
                goto bisect;
            } else if (ok1 && ok2) {
                d = (std::fabs(d1) < std::fabs(d2)) ? d1 : d2;
            } else if (ok1) {
                d = d1;
            } else {
                d = d2;
            }

            if (std::fabs(d) > std::fabs(0.5 * olde)) goto bisect;

            {
                double u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = std::copysign(tol1, xm - x);
                goto evaluate;
            }
        }

    bisect:
        e = (dx >= 0.0) ? a - x : b - x;
        d = 0.5 * e;

    evaluate:
        double u, fu;
        if (std::fabs(d) >= tol1) {
            u  = x + d;
            fu = f(u);
        } else {
            u  = x + std::copysign(tol1, d);
            fu = f(u);
            if (fu > fx) {
                xmin = x;
                return fx;
            }
        }

        double du = df(u);
        if (fu <= fx) {
            if (u >= x) a = x; else b = x;
            v = w;  fv = fw;  dv = dw;
            w = x;  fw = fx;  dw = dx;
            x = u;  fx = fu;  dx = du;
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w;  fv = fw;  dv = dw;
                w = u;  fw = fu;  dw = du;
            } else if (fu <= fv || v == x || v == w) {
                v = u;  fv = fu;  dv = du;
            }
        }
    }
    std::fprintf(stderr, "dbrent: exceeded maximum iterations\n");
    xmin = x;
    return fx;
}
