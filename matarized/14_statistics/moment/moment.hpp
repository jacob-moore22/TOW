#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

// Compute mean, average deviation, standard deviation, variance,
// skewness, and kurtosis of data(1..n).
inline void moment(DFMatrixKokkos<double>& data, int n,
                   double& ave, double& adev, double& sdev,
                   double& var, double& skew, double& curt)
{
    if (n <= 1) {
        std::printf("moment: N must be at least 2\n");
        return;
    }

    double s = 0.0;
    DO_REDUCE_SUM(j, 1, n, s, {
        s += data(j);
    });
    ave = s / n;

    adev = 0.0;
    var  = 0.0;
    skew = 0.0;
    curt = 0.0;

    double sum_adev = 0.0, sum_var = 0.0, sum_skew = 0.0, sum_curt = 0.0;
    for (int j = 1; j <= n; j++) {
        double d = data.host(j) - ave;
        sum_adev += fabs(d);
        double p = d * d;
        sum_var += p;
        p *= d;
        sum_skew += p;
        p *= d;
        sum_curt += p;
    }
    adev = sum_adev / n;
    var  = sum_var / (n - 1);
    sdev = sqrt(var);

    if (var != 0.0) {
        skew = sum_skew / (n * sdev * sdev * sdev);
        curt = sum_curt / (n * var * var) - 3.0;
    } else {
        std::printf("moment: no skew or kurtosis when zero variance\n");
    }
}
