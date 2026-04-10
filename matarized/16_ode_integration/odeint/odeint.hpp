// ODE integrator driver with adaptive stepping (Numerical Recipes ODEINT).
// Integrates y[1..nvar] from x1 to x2 using a user-supplied stepper
// (rkqc or bsstep). COMMON /PATH/ converted to OdeIntPath struct for
// optional storage of intermediate results.

#pragma once
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <matar.h>

using namespace mtr;

struct OdeIntPath {
    int kmax    = 0;
    double dxsav = 0.0;
    int kount   = 0;
    DFMatrixKokkos<double> xp;
    DFMatrixKokkos<double> yp;
};

template<typename Derivs, typename Stepper>
inline void odeint(DFMatrixKokkos<double>& ystart, int nvar,
                   double x1, double x2, double eps, double h1, double hmin,
                   int& nok, int& nbad,
                   Derivs derivs, Stepper stepper,
                   OdeIntPath& path)
{
    constexpr int    MAXSTP = 10000;
    constexpr double TINY   = 1.0e-30;

    DFMatrixKokkos<double> yscal(nvar), y(nvar), dydx(nvar);

    double x = x1;
    double h = std::copysign(h1, x2 - x1);
    nok  = 0;
    nbad = 0;
    path.kount = 0;

    for (int i = 1; i <= nvar; i++)
        y.host(i) = ystart.host(i);

    double xsav = x - path.dxsav * 2.0;

    for (int nstp = 1; nstp <= MAXSTP; nstp++) {
        derivs(x, y, dydx);
        for (int i = 1; i <= nvar; i++)
            yscal.host(i) = std::fabs(y.host(i))
                          + std::fabs(h * dydx.host(i)) + TINY;

        if (path.kmax > 0) {
            if (std::fabs(x - xsav) > std::fabs(path.dxsav)) {
                if (path.kount < path.kmax - 1) {
                    path.kount++;
                    path.xp.host(path.kount) = x;
                    for (int i = 1; i <= nvar; i++)
                        path.yp.host(i, path.kount) = y.host(i);
                    xsav = x;
                }
            }
        }

        if ((x + h - x2) * (x + h - x1) > 0.0)
            h = x2 - x;

        double hdid, hnext;
        stepper(y, dydx, nvar, x, h, eps, yscal, hdid, hnext, derivs);

        if (hdid == h)
            nok++;
        else
            nbad++;

        if ((x - x2) * (x2 - x1) >= 0.0) {
            for (int i = 1; i <= nvar; i++)
                ystart.host(i) = y.host(i);
            if (path.kmax != 0) {
                path.kount++;
                path.xp.host(path.kount) = x;
                for (int i = 1; i <= nvar; i++)
                    path.yp.host(i, path.kount) = y.host(i);
            }
            return;
        }

        if (std::fabs(hnext) < hmin) {
            std::fprintf(stderr, "odeint: stepsize smaller than minimum\n");
            return;
        }
        h = hnext;
    }
    std::fprintf(stderr, "odeint: too many steps\n");
}
