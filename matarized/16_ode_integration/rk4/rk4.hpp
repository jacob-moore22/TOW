// Classic 4th-order Runge-Kutta single step (Numerical Recipes RK4).
// Advances y[1..n] from x by step h, given derivatives dydx at x.
// y and yout may alias (in-place update is safe).

#pragma once
#include <matar.h>

using namespace mtr;

template<typename Derivs>
inline void rk4(DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& dydx,
                int n, double x, double h, DFMatrixKokkos<double>& yout,
                Derivs derivs)
{
    DFMatrixKokkos<double> yt(n), dyt(n), dym(n);

    double hh = h * 0.5;
    double h6 = h / 6.0;
    double xh = x + hh;

    for (int i = 1; i <= n; i++)
        yt.host(i) = y.host(i) + hh * dydx.host(i);

    derivs(xh, yt, dyt);

    for (int i = 1; i <= n; i++)
        yt.host(i) = y.host(i) + hh * dyt.host(i);

    derivs(xh, yt, dym);

    for (int i = 1; i <= n; i++) {
        yt.host(i) = y.host(i) + h * dym.host(i);
        dym.host(i) = dyt.host(i) + dym.host(i);
    }

    derivs(x + h, yt, dyt);

    for (int i = 1; i <= n; i++)
        yout.host(i) = y.host(i) + h6 * (dydx.host(i) + dyt.host(i) + 2.0 * dym.host(i));
}
