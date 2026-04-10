#pragma once
#include <cmath>
#include <algorithm>
#include <matar.h>
#include "gammq.hpp"

using namespace mtr;

// Chi-square test on a contingency table nn(1..ni, 1..nj).
// Returns chisq, df, prob, Cramer's V (cramrv), and contingency coefficient (ccc).
inline void cntab1(DFMatrixKokkos<double>& nn, int ni, int nj,
                   double& chisq, double& df, double& prob,
                   double& cramrv, double& ccc)
{
    constexpr double TINY = 1.0e-30;

    DFMatrixKokkos<double> sumi(ni), sumj(nj);

    double sum = 0.0;
    int nni = ni;
    int nnj = nj;

    for (int i = 1; i <= ni; i++) {
        sumi.host(i) = 0.0;
        for (int j = 1; j <= nj; j++) {
            sumi.host(i) += nn.host(i, j);
            sum += nn.host(i, j);
        }
        if (sumi.host(i) == 0.0) nni--;
    }

    for (int j = 1; j <= nj; j++) {
        sumj.host(j) = 0.0;
        for (int i = 1; i <= ni; i++)
            sumj.host(j) += nn.host(i, j);
        if (sumj.host(j) == 0.0) nnj--;
    }

    df = nni * nnj - nni - nnj + 1;
    chisq = 0.0;
    for (int i = 1; i <= ni; i++) {
        for (int j = 1; j <= nj; j++) {
            double expctd = sumj.host(j) * sumi.host(i) / sum;
            chisq += (nn.host(i, j) - expctd) * (nn.host(i, j) - expctd) / (expctd + TINY);
        }
    }

    prob = gammq(0.5 * df, 0.5 * chisq);
    cramrv = sqrt(chisq / (sum * std::min(nni - 1, nnj - 1)));
    ccc = sqrt(chisq / (chisq + sum));
}
