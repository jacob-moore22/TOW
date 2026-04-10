// Fixed-step RK4 driver (Numerical Recipes RKDUMB).
// Integrates y[1..nvar] from x1 to x2 using nstep equal RK4 steps.
// COMMON /PATH/ converted to output arrays xx[1..nstep+1] and y_out(nvar,nstep+1).

#pragma once
#include <cstdio>
#include <matar.h>
#include "../rk4/rk4.hpp"

using namespace mtr;

template<typename Derivs>
inline void rkdumb(DFMatrixKokkos<double>& vstart, int nvar, double x1, double x2,
                   int nstep, Derivs derivs,
                   DFMatrixKokkos<double>& xx, DFMatrixKokkos<double>& y_out)
{
    DFMatrixKokkos<double> v(nvar), dv(nvar);

    for (int i = 1; i <= nvar; i++) {
        v.host(i) = vstart.host(i);
        y_out.host(i, 1) = v.host(i);
    }

    xx.host(1) = x1;
    double x = x1;
    double h = (x2 - x1) / nstep;

    for (int k = 1; k <= nstep; k++) {
        derivs(x, v, dv);
        rk4(v, dv, nvar, x, h, v, derivs);
        if (x + h == x) {
            std::fprintf(stderr, "rkdumb: stepsize not significant\n");
            return;
        }
        x += h;
        xx.host(k + 1) = x;
        for (int i = 1; i <= nvar; i++)
            y_out.host(i, k + 1) = v.host(i);
    }
}
