#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

// Golden section search for a minimum of f(x) in the bracket [ax, cx] with
// interior point bx. Returns the minimum function value; xmin is set to the
// location of the minimum.
template <typename Func>
inline double golden(double ax, double bx, double cx, Func f,
                     double tol, double& xmin)
{
    constexpr double R = 0.61803399;
    constexpr double C = 0.38196602;

    double x0 = ax;
    double x3 = cx;
    double x1, x2;

    if (std::fabs(cx - bx) > std::fabs(bx - ax)) {
        x1 = bx;
        x2 = bx + C * (cx - bx);
    } else {
        x2 = bx;
        x1 = bx - C * (bx - ax);
    }

    double f1 = f(x1);
    double f2 = f(x2);

    while (std::fabs(x3 - x0) > tol * (std::fabs(x1) + std::fabs(x2))) {
        if (f2 < f1) {
            x0 = x1;
            x1 = x2;
            x2 = R * x1 + C * x3;
            f1 = f2;
            f2 = f(x2);
        } else {
            x3 = x2;
            x2 = x1;
            x1 = R * x2 + C * x0;
            f2 = f1;
            f1 = f(x1);
        }
    }

    if (f1 < f2) {
        xmin = x1;
        return f1;
    } else {
        xmin = x2;
        return f2;
    }
}
