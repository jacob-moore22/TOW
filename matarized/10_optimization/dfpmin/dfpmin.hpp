#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "../linmin/linmin.hpp"

using namespace mtr;

// BFGS quasi-Newton minimization.
//
// p:     starting point (length n), replaced by the minimum on return.
// n:     number of dimensions.
// ftol:  fractional convergence tolerance on the function value.
// iter:  on return, number of iterations performed.
// fret:  on return, minimum function value.
// func:  objective  f(const double* x) -> double.
// dfunc: gradient  dfunc(const double* x, double* grad).
template <typename Func, typename DFunc>
inline void dfpmin(double* p, int n, double ftol, int& iter, double& fret,
                   Func func, DFunc dfunc)
{
    constexpr int    NMAX  = 50;
    constexpr int    ITMAX = 200;
    constexpr double EPS   = 1.0e-10;

    double hessin[NMAX * NMAX], xi[NMAX], g[NMAX], dg[NMAX], hdg[NMAX];

    double fp = func(p);
    dfunc(p, g);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            hessin[i * n + j] = 0.0;
        hessin[i * n + i] = 1.0;
        xi[i] = -g[i];
    }

    for (int its = 0; its < ITMAX; its++) {
        iter = its + 1;
        linmin(p, xi, n, fret, func);

        if (2.0 * std::fabs(fret - fp) <= ftol * (std::fabs(fret) + std::fabs(fp) + EPS))
            return;

        fp = fret;
        for (int i = 0; i < n; i++) dg[i] = g[i];

        fret = func(p);
        dfunc(p, g);
        for (int i = 0; i < n; i++) dg[i] = g[i] - dg[i];

        for (int i = 0; i < n; i++) {
            hdg[i] = 0.0;
            for (int j = 0; j < n; j++)
                hdg[i] += hessin[i * n + j] * dg[j];
        }

        double fac = 0.0, fae = 0.0;
        for (int i = 0; i < n; i++) {
            fac += dg[i] * xi[i];
            fae += dg[i] * hdg[i];
        }
        fac = 1.0 / fac;
        double fad = 1.0 / fae;

        for (int i = 0; i < n; i++)
            dg[i] = fac * xi[i] - fad * hdg[i];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                hessin[i * n + j] += fac * xi[i] * xi[j]
                                   - fad * hdg[i] * hdg[j]
                                   + fae * dg[i] * dg[j];
            }
        }

        for (int i = 0; i < n; i++) {
            xi[i] = 0.0;
            for (int j = 0; j < n; j++)
                xi[i] -= hessin[i * n + j] * g[j];
        }
    }
    std::fprintf(stderr, "dfpmin: exceeded maximum iterations\n");
}
