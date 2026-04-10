#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

// Compute mean and variance of data(1..n).
inline void avevar(DFMatrixKokkos<double>& data, int n,
                   double& ave, double& var)
{
    double s = 0.0;
    DO_REDUCE_SUM(j, 1, n, s, {
        s += data(j);
    });
    ave = s / n;

    double sv = 0.0;
    for (int j = 1; j <= n; j++) {
        double d = data.host(j) - ave;
        sv += d * d;
    }
    var = sv / (n - 1);
}
