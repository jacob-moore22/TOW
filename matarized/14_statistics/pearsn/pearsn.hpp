#pragma once
#include <cmath>
#include <matar.h>
#include "betai.hpp"

using namespace mtr;

// Pearson's correlation coefficient r, its significance prob,
// and Fisher's z for arrays x(1..n) and y(1..n).
inline void pearsn(DFMatrixKokkos<double>& x, DFMatrixKokkos<double>& y,
                   int n, double& r, double& prob, double& z)
{
    constexpr double TINY = 1.0e-20;

    double ax = 0.0, ay = 0.0;
    for (int j = 1; j <= n; j++) {
        ax += x.host(j);
        ay += y.host(j);
    }
    ax /= n;
    ay /= n;

    double sxx = 0.0, syy = 0.0, sxy = 0.0;
    for (int j = 1; j <= n; j++) {
        double xt = x.host(j) - ax;
        double yt = y.host(j) - ay;
        sxx += xt * xt;
        syy += yt * yt;
        sxy += xt * yt;
    }

    r = sxy / sqrt(sxx * syy);
    z = 0.5 * log(((1.0 + r) + TINY) / ((1.0 - r) + TINY));
    double df = n - 2;
    double t = r * sqrt(df / (((1.0 - r) + TINY) * ((1.0 + r) + TINY)));
    prob = betai(0.5 * df, 0.5, df / (df + t * t));
}
