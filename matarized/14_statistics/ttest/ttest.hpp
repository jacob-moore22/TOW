#pragma once
#include <cmath>
#include <matar.h>
#include "avevar.hpp"
#include "betai.hpp"

using namespace mtr;

// Student's t-test for equal variances.
// Given data1(1..n1) and data2(1..n2), returns the t statistic and
// significance level prob (small prob = significant difference).
inline void ttest(DFMatrixKokkos<double>& data1, int n1,
                  DFMatrixKokkos<double>& data2, int n2,
                  double& t, double& prob)
{
    double ave1, var1, ave2, var2;
    avevar(data1, n1, ave1, var1);
    avevar(data2, n2, ave2, var2);

    double df = n1 + n2 - 2;
    double var = ((n1 - 1) * var1 + (n2 - 1) * var2) / df;
    t = (ave1 - ave2) / sqrt(var * (1.0 / n1 + 1.0 / n2));
    prob = betai(0.5 * df, 0.5, df / (df + t * t));
}
