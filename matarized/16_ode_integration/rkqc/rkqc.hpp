// RK4 with adaptive step control via Richardson extrapolation
// (Numerical Recipes RKQC).
// Takes a single quality-controlled Runge-Kutta step.
// Compares two half-steps with one full step to estimate error.

#pragma once
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <matar.h>
#include "../rk4/rk4.hpp"

using namespace mtr;

template<typename Derivs>
inline void rkqc(DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& dydx,
                 int n, double& x, double htry, double eps,
                 DFMatrixKokkos<double>& yscal, double& hdid, double& hnext,
                 Derivs derivs)
{
    constexpr double FCOR   = 1.0 / 15.0;
    constexpr double SAFETY = 0.9;
    constexpr double ERRCON = 6.0e-4;

    double pgrow  = -0.20;
    double pshrnk = -0.25;
    double xsav   = x;

    DFMatrixKokkos<double> ytemp(n), ysav(n), dysav(n);

    for (int i = 1; i <= n; i++) {
        ysav.host(i)  = y.host(i);
        dysav.host(i) = dydx.host(i);
    }

    double h = htry;
    for (;;) {
        double hh = 0.5 * h;

        rk4(ysav, dysav, n, xsav, hh, ytemp, derivs);
        x = xsav + hh;
        derivs(x, ytemp, dydx);
        rk4(ytemp, dydx, n, x, hh, y, derivs);

        x = xsav + h;
        if (x == xsav) {
            std::fprintf(stderr, "rkqc: stepsize not significant\n");
            return;
        }

        rk4(ysav, dysav, n, xsav, h, ytemp, derivs);

        double errmax = 0.0;
        for (int i = 1; i <= n; i++) {
            ytemp.host(i) = y.host(i) - ytemp.host(i);
            errmax = std::max(errmax,
                              std::fabs(ytemp.host(i) / yscal.host(i)));
        }
        errmax /= eps;

        if (errmax > 1.0) {
            h = SAFETY * h * std::pow(errmax, pshrnk);
        } else {
            hdid = h;
            hnext = (errmax > ERRCON)
                        ? SAFETY * h * std::pow(errmax, pgrow)
                        : 4.0 * h;
            break;
        }
    }

    for (int i = 1; i <= n; i++)
        y.host(i) += ytemp.host(i) * FCOR;
}
