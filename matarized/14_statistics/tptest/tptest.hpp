#pragma once
#include <cmath>
#include <matar.h>
#include "avevar.hpp"
#include "betai.hpp"

using namespace mtr;

// Paired-sample t-test for data1(1..n) and data2(1..n).
inline void tptest(DFMatrixKokkos<double>& data1,
                   DFMatrixKokkos<double>& data2, int n,
                   double& t, double& prob)
{
    double ave1, var1, ave2, var2;
    avevar(data1, n, ave1, var1);
    avevar(data2, n, ave2, var2);

    double cov = 0.0;
    for (int j = 1; j <= n; j++)
        cov += (data1.host(j) - ave1) * (data2.host(j) - ave2);

    double df = n - 1;
    cov /= df;
    double sd = sqrt((var1 + var2 - 2.0 * cov) / n);
    t = (ave1 - ave2) / sd;
    prob = betai(0.5 * df, 0.5, df / (df + t * t));
}
