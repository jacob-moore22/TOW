// Bulirsch-Stoer adaptive step (Numerical Recipes BSSTEP).
// Takes a single quality-controlled step using modified midpoint + rational
// extrapolation. Depends on: mmid, rzextr.

#pragma once
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <matar.h>
#include "../mmid/mmid.hpp"
#include "../rzextr/rzextr.hpp"

using namespace mtr;

template<typename Derivs>
inline void bsstep(DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& dydx,
                   int nv, double& x, double htry, double eps,
                   DFMatrixKokkos<double>& yscal, double& hdid, double& hnext,
                   Derivs derivs)
{
    constexpr int    IMAX   = 11;
    constexpr int    NUSE   = 7;
    constexpr double SHRINK = 0.95;
    constexpr double GROW   = 1.2;

    static constexpr int nseq[] = {0, 2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96};

    DFMatrixKokkos<double> yerr(nv), ysav(nv), dysav(nv), yseq(nv);

    double h    = htry;
    double xsav = x;

    for (int i = 1; i <= nv; i++) {
        ysav.host(i)  = y.host(i);
        dysav.host(i) = dydx.host(i);
    }

    for (;;) {
        for (int ii = 1; ii <= IMAX; ii++) {
            mmid(ysav, dysav, nv, xsav, h, nseq[ii], yseq, derivs);
            double xest = (h / nseq[ii]) * (h / nseq[ii]);
            rzextr(ii, xest, yseq, y, yerr, nv, NUSE);

            double errmax = 0.0;
            for (int j = 1; j <= nv; j++)
                errmax = std::max(errmax,
                                  std::fabs(yerr.host(j) / yscal.host(j)));
            errmax /= eps;

            if (errmax < 1.0) {
                x    = xsav + h;
                hdid = h;
                if (ii == NUSE)
                    hnext = h * SHRINK;
                else if (ii == NUSE - 1)
                    hnext = h * GROW;
                else
                    hnext = (h * nseq[NUSE - 1]) / nseq[ii];
                return;
            }
        }

        h *= 0.25 / (1 << ((IMAX - NUSE) / 2));
        if (x + h == x) {
            std::fprintf(stderr, "bsstep: step size underflow\n");
            return;
        }
    }
}
