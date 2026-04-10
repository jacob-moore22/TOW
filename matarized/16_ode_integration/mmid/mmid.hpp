// Modified midpoint method (Numerical Recipes MMID).
// Advances y[1..nvar] from xs by total step htot using nstep substeps.
// Used as the base integrator for the Bulirsch-Stoer method.

#pragma once
#include <matar.h>

using namespace mtr;

template<typename Derivs>
inline void mmid(DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& dydx,
                 int nvar, double xs, double htot, int nstep,
                 DFMatrixKokkos<double>& yout, Derivs derivs)
{
    DFMatrixKokkos<double> ym(nvar), yn(nvar);

    double h = htot / nstep;
    for (int i = 1; i <= nvar; i++) {
        ym.host(i) = y.host(i);
        yn.host(i) = y.host(i) + h * dydx.host(i);
    }

    double x = xs + h;
    derivs(x, yn, yout);

    double h2 = 2.0 * h;
    for (int ns = 2; ns <= nstep; ns++) {
        for (int i = 1; i <= nvar; i++) {
            double swap = ym.host(i) + h2 * yout.host(i);
            ym.host(i) = yn.host(i);
            yn.host(i) = swap;
        }
        x += h;
        derivs(x, yn, yout);
    }

    for (int i = 1; i <= nvar; i++)
        yout.host(i) = 0.5 * (ym.host(i) + yn.host(i) + h * yout.host(i));
}
