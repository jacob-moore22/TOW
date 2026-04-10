#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "../linmin/linmin.hpp"

using namespace mtr;

// Fletcher-Reeves-Polak-Ribiere conjugate gradient minimization.
//
// p:     starting point (length n), replaced by the minimum on return.
// n:     number of dimensions.
// ftol:  fractional convergence tolerance on the function value.
// iter:  on return, number of iterations performed.
// fret:  on return, minimum function value.
// func:  objective  f(const double* x) -> double.
// dfunc: gradient  dfunc(const double* x, double* grad).
template <typename Func, typename DFunc>
inline void frprmn(double* p, int n, double ftol, int& iter, double& fret,
                   Func func, DFunc dfunc)
{
    constexpr int    NMAX  = 50;
    constexpr int    ITMAX = 200;
    constexpr double EPS   = 1.0e-10;

    double g[NMAX], h[NMAX], xi[NMAX];

    double fp = func(p);
    dfunc(p, xi);
    for (int j = 0; j < n; j++) {
        g[j]  = -xi[j];
        h[j]  =  g[j];
        xi[j] =  h[j];
    }

    for (int its = 0; its < ITMAX; its++) {
        iter = its + 1;
        linmin(p, xi, n, fret, func);

        if (2.0 * std::fabs(fret - fp) <= ftol * (std::fabs(fret) + std::fabs(fp) + EPS))
            return;

        fp = func(p);
        dfunc(p, xi);

        double gg  = 0.0;
        double dgg = 0.0;
        for (int j = 0; j < n; j++) {
            gg  += g[j] * g[j];
            // Polak-Ribiere formula (comment out next line and uncomment
            // the one after for Fletcher-Reeves)
            dgg += (xi[j] + g[j]) * xi[j];
            // dgg += xi[j] * xi[j];
        }
        if (gg == 0.0) return;

        double gam = dgg / gg;
        for (int j = 0; j < n; j++) {
            g[j]  = -xi[j];
            h[j]  = g[j] + gam * h[j];
            xi[j] = h[j];
        }
    }
    std::fprintf(stderr, "frprmn: exceeded maximum iterations\n");
}
