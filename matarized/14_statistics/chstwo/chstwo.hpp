#pragma once
#include <cmath>
#include <matar.h>
#include "gammq.hpp"

using namespace mtr;

// Chi-square two-sample test. Compares bins1(1..nbins) vs bins2(1..nbins).
// knstrn = number of constraints.
inline void chstwo(DFMatrixKokkos<double>& bins1, DFMatrixKokkos<double>& bins2,
                   int nbins, int knstrn,
                   double& df, double& chsq, double& prob)
{
    df = nbins - 1 - knstrn;
    chsq = 0.0;
    for (int j = 1; j <= nbins; j++) {
        if (bins1.host(j) == 0.0 && bins2.host(j) == 0.0) {
            df -= 1.0;
        } else {
            double diff = bins1.host(j) - bins2.host(j);
            chsq += diff * diff / (bins1.host(j) + bins2.host(j));
        }
    }
    prob = gammq(0.5 * df, 0.5 * chsq);
}
