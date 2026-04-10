#pragma once
#include <cmath>
#include <matar.h>
#include "sort2.hpp"
#include "crank.hpp"
#include "betai.hpp"
#include "erfcc.hpp"

using namespace mtr;

// Spearman rank-order correlation coefficient.
// Given data1(1..n) and data2(1..n), workspace wksp1(1..n) and wksp2(1..n).
// Returns: d (sum of squared rank differences), zd, probd (significance via normal approx),
// rs (Spearman's r_s), probrs (significance via t-test with betai).
inline void spear(DFMatrixKokkos<double>& data1, DFMatrixKokkos<double>& data2,
                  int n,
                  DFMatrixKokkos<double>& wksp1, DFMatrixKokkos<double>& wksp2,
                  double& d, double& zd, double& probd,
                  double& rs, double& probrs)
{
    for (int j = 1; j <= n; j++) {
        wksp1.host(j) = data1.host(j);
        wksp2.host(j) = data2.host(j);
    }

    sort2(n, wksp1, wksp2);
    double sf;
    crank(n, wksp1, sf);

    sort2(n, wksp2, wksp1);
    double sg;
    crank(n, wksp2, sg);

    d = 0.0;
    for (int j = 1; j <= n; j++) {
        double diff = wksp1.host(j) - wksp2.host(j);
        d += diff * diff;
    }

    double en   = n;
    double en3n = en * en * en - en;
    double aved = en3n / 6.0 - (sf + sg) / 12.0;
    double fac  = (1.0 - sf / en3n) * (1.0 - sg / en3n);
    double vard = ((en - 1.0) * en * en * (en + 1.0) * (en + 1.0) / 36.0) * fac;

    zd    = (d - aved) / sqrt(vard);
    probd = erfcc(fabs(zd) / 1.4142136);

    rs = (1.0 - (6.0 / en3n) * (d + 0.5 * (sf + sg))) / fac;
    double t  = rs * sqrt((en - 2.0) / ((1.0 + rs) * (1.0 - rs)));
    double df = en - 2.0;
    probrs = betai(0.5 * df, 0.5, df / (df + t * t));
}
