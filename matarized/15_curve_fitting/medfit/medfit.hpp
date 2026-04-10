// Median-based robust line fitting (Numerical Recipes MEDFIT).
//
// Fits y = a + b*x by minimizing the sum of absolute deviations
// (L1 norm), which is resistant to outliers. Uses bisection on
// rofunc to find the optimal slope, then computes the intercept
// as the median of residuals.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "rofunc.hpp"

using namespace mtr;

inline void medfit(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                   int ndata, double& a, double& b, double& abdev)
{
    double sx = 0.0, sy = 0.0, sxy = 0.0, sxx = 0.0;

    DFMatrixKokkos<double> xt(ndata), yt(ndata), arr(ndata);

    for (int j = 1; j <= ndata; j++) {
        xt.host(j) = x.host(j);
        yt.host(j) = y.host(j);
        sx  += x.host(j);
        sy  += y.host(j);
        sxy += x.host(j) * y.host(j);
        sxx += x.host(j) * x.host(j);
    }

    double del = static_cast<double>(ndata) * sxx - sx * sx;
    double aa  = (sxx * sy - sx * sxy) / del;
    double bb  = (static_cast<double>(ndata) * sxy - sx * sy) / del;

    double chisq = 0.0;
    for (int j = 1; j <= ndata; j++) {
        double r = y.host(j) - (aa + bb * x.host(j));
        chisq += r * r;
    }

    double sigb = std::sqrt(chisq / del);
    double abdevt;

    double b1 = bb;
    double f1 = rofunc(b1, ndata, xt, yt, arr, aa, abdevt);
    double b2 = bb + std::copysign(3.0 * sigb, f1);
    double f2 = rofunc(b2, ndata, xt, yt, arr, aa, abdevt);

    while (f1 * f2 > 0.0) {
        bb = 2.0 * b2 - b1;
        b1 = b2;
        f1 = f2;
        b2 = bb;
        f2 = rofunc(b2, ndata, xt, yt, arr, aa, abdevt);
    }

    sigb *= 0.01;
    while (std::fabs(b2 - b1) > sigb) {
        bb = 0.5 * (b1 + b2);
        if (bb == b1 || bb == b2) break;
        double f = rofunc(bb, ndata, xt, yt, arr, aa, abdevt);
        if (f * f1 >= 0.0) {
            f1 = f;
            b1 = bb;
        } else {
            f2 = f;
            b2 = bb;
        }
    }

    a = aa;
    b = bb;
    abdev = abdevt / ndata;
}
