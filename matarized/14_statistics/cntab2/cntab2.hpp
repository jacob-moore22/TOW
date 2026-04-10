#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

// Entropy-based measures of association for contingency table nn(1..ni, 1..nj).
// Returns H (joint entropy), HX, HY (marginal entropies),
// HYGX, HXGY (conditional entropies), UYGX, UXGY, UXY (uncertainty coefficients).
inline void cntab2(DFMatrixKokkos<double>& nn, int ni, int nj,
                   double& h, double& hx, double& hy,
                   double& hygx, double& hxgy,
                   double& uygx, double& uxgy, double& uxy)
{
    constexpr double TINY = 1.0e-30;

    DFMatrixKokkos<double> sumi(ni), sumj(nj);

    double sum = 0.0;
    for (int i = 1; i <= ni; i++) {
        sumi.host(i) = 0.0;
        for (int j = 1; j <= nj; j++) {
            sumi.host(i) += nn.host(i, j);
            sum += nn.host(i, j);
        }
    }
    for (int j = 1; j <= nj; j++) {
        sumj.host(j) = 0.0;
        for (int i = 1; i <= ni; i++)
            sumj.host(j) += nn.host(i, j);
    }

    hx = 0.0;
    for (int i = 1; i <= ni; i++) {
        if (sumi.host(i) != 0.0) {
            double p = sumi.host(i) / sum;
            hx -= p * log(p);
        }
    }

    hy = 0.0;
    for (int j = 1; j <= nj; j++) {
        if (sumj.host(j) != 0.0) {
            double p = sumj.host(j) / sum;
            hy -= p * log(p);
        }
    }

    h = 0.0;
    for (int i = 1; i <= ni; i++) {
        for (int j = 1; j <= nj; j++) {
            if (nn.host(i, j) != 0.0) {
                double p = nn.host(i, j) / sum;
                h -= p * log(p);
            }
        }
    }

    hygx = h - hx;
    hxgy = h - hy;
    uygx = (hy - hygx) / (hy + TINY);
    uxgy = (hx - hxgy) / (hx + TINY);
    uxy  = 2.0 * (hx + hy - h) / (hx + hy + TINY);
}
