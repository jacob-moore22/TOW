#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Nelder-Mead downhill simplex method. Minimizes an n-dimensional function.
//
// p:     simplex vertices stored row-major in a flat array of size (ndim+1)*ndim.
//        p[i*ndim + j] is the j-th coordinate of vertex i (0-based).
// y:     function values at each vertex, length ndim+1.
// ndim:  number of dimensions.
// ftol:  fractional convergence tolerance.
// funk:  objective function taking (const double* x) -> double.
// iter:  on return, number of function evaluations used.
template <typename Func>
inline void amoeba(double* p, double* y, int ndim, double ftol,
                   Func funk, int& iter)
{
    constexpr double ALPHA = 1.0;
    constexpr double BETA  = 0.5;
    constexpr double GAMMA = 2.0;
    constexpr int    ITMAX = 5000;
    constexpr int    NMAX  = 50;

    double pr[NMAX], prr[NMAX], pbar[NMAX];

    int mpts = ndim + 1;
    iter = 0;

    for (;;) {
        int ilo = 0;
        int ihi, inhi;
        if (y[0] > y[1]) { ihi = 0; inhi = 1; }
        else              { ihi = 1; inhi = 0; }

        for (int i = 0; i < mpts; i++) {
            if (y[i] < y[ilo]) ilo = i;
            if (y[i] > y[ihi]) {
                inhi = ihi;
                ihi  = i;
            } else if (y[i] > y[inhi]) {
                if (i != ihi) inhi = i;
            }
        }

        double rtol = 2.0 * std::fabs(y[ihi] - y[ilo]) /
                      (std::fabs(y[ihi]) + std::fabs(y[ilo]) + 1.0e-30);
        if (rtol < ftol) return;
        if (iter >= ITMAX) {
            std::fprintf(stderr, "amoeba: exceeded maximum iterations\n");
            return;
        }
        iter++;

        for (int j = 0; j < ndim; j++) pbar[j] = 0.0;

        for (int i = 0; i < mpts; i++) {
            if (i != ihi) {
                for (int j = 0; j < ndim; j++)
                    pbar[j] += p[i * ndim + j];
            }
        }

        for (int j = 0; j < ndim; j++) {
            pbar[j] /= ndim;
            pr[j] = (1.0 + ALPHA) * pbar[j] - ALPHA * p[ihi * ndim + j];
        }

        double ypr = funk(pr);

        if (ypr <= y[ilo]) {
            for (int j = 0; j < ndim; j++)
                prr[j] = GAMMA * pr[j] + (1.0 - GAMMA) * pbar[j];
            double yprr = funk(prr);
            if (yprr < y[ilo]) {
                for (int j = 0; j < ndim; j++) p[ihi * ndim + j] = prr[j];
                y[ihi] = yprr;
            } else {
                for (int j = 0; j < ndim; j++) p[ihi * ndim + j] = pr[j];
                y[ihi] = ypr;
            }
        } else if (ypr >= y[inhi]) {
            if (ypr < y[ihi]) {
                for (int j = 0; j < ndim; j++) p[ihi * ndim + j] = pr[j];
                y[ihi] = ypr;
            }
            for (int j = 0; j < ndim; j++)
                prr[j] = BETA * p[ihi * ndim + j] + (1.0 - BETA) * pbar[j];
            double yprr = funk(prr);
            if (yprr < y[ihi]) {
                for (int j = 0; j < ndim; j++) p[ihi * ndim + j] = prr[j];
                y[ihi] = yprr;
            } else {
                for (int i = 0; i < mpts; i++) {
                    if (i != ilo) {
                        for (int j = 0; j < ndim; j++) {
                            pr[j] = 0.5 * (p[i * ndim + j] + p[ilo * ndim + j]);
                            p[i * ndim + j] = pr[j];
                        }
                        y[i] = funk(pr);
                    }
                }
            }
        } else {
            for (int j = 0; j < ndim; j++) p[ihi * ndim + j] = pr[j];
            y[ihi] = ypr;
        }
    }
}
