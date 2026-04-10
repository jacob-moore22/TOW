#pragma once
#include <cstdio>
#include <cmath>
#include <matar.h>
#include "gammq.hpp"

using namespace mtr;

// Chi-square one-sample test. Compares observed bins(1..nbins) to
// expected ebins(1..nbins). knstrn = number of fitted parameters.
inline void chsone(DFMatrixKokkos<double>& bins, DFMatrixKokkos<double>& ebins,
                   int nbins, int knstrn,
                   double& df, double& chsq, double& prob)
{
    df = nbins - 1 - knstrn;
    chsq = 0.0;
    for (int j = 1; j <= nbins; j++) {
        if (ebins.host(j) <= 0.0) {
            std::printf("chsone: bad expected number in bin %d\n", j);
            return;
        }
        double diff = bins.host(j) - ebins.host(j);
        chsq += diff * diff / ebins.host(j);
    }
    prob = gammq(0.5 * df, 0.5 * chsq);
}
