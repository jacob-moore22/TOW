// Linear least squares fit y = a + bx (Numerical Recipes FIT).
//
// Fits data (x, y) with optional per-point weights (sig).
// Returns: intercept a, slope b, their uncertainties siga/sigb,
// chi-squared, and goodness-of-fit probability q.
//
// mwt == 0: unweighted (sig ignored, q set to 1).
// mwt != 0: weighted by 1/sig(i)^2.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "gammq.hpp"

using namespace mtr;

inline void fit(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                int ndata, DFMatrixKokkos<double>& sig, int mwt,
                double& a, double& b, double& siga, double& sigb,
                double& chi2, double& q)
{
    double sx = 0.0, sy = 0.0, st2 = 0.0, ss;
    b = 0.0;

    if (mwt != 0) {
        ss = 0.0;
        for (int i = 1; i <= ndata; i++) {
            double wt = 1.0 / (sig.host(i) * sig.host(i));
            ss += wt;
            sx += x.host(i) * wt;
            sy += y.host(i) * wt;
        }
    } else {
        for (int i = 1; i <= ndata; i++) {
            sx += x.host(i);
            sy += y.host(i);
        }
        ss = static_cast<double>(ndata);
    }

    double sxoss = sx / ss;

    if (mwt != 0) {
        for (int i = 1; i <= ndata; i++) {
            double t = (x.host(i) - sxoss) / sig.host(i);
            st2 += t * t;
            b += t * y.host(i) / sig.host(i);
        }
    } else {
        for (int i = 1; i <= ndata; i++) {
            double t = x.host(i) - sxoss;
            st2 += t * t;
            b += t * y.host(i);
        }
    }

    b /= st2;
    a = (sy - sx * b) / ss;
    siga = std::sqrt((1.0 + sx * sx / (ss * st2)) / ss);
    sigb = std::sqrt(1.0 / st2);

    chi2 = 0.0;
    if (mwt == 0) {
        for (int i = 1; i <= ndata; i++) {
            double r = y.host(i) - a - b * x.host(i);
            chi2 += r * r;
        }
        q = 1.0;
        double sigdat = std::sqrt(chi2 / (ndata - 2));
        siga *= sigdat;
        sigb *= sigdat;
    } else {
        for (int i = 1; i <= ndata; i++) {
            double r = (y.host(i) - a - b * x.host(i)) / sig.host(i);
            chi2 += r * r;
        }
        q = gammq(0.5 * (ndata - 2), 0.5 * chi2);
    }
}
