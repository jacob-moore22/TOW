#pragma once
#include <cmath>
#include <matar.h>
#include "avevar.hpp"
#include "betai.hpp"

using namespace mtr;

// Student's t-test for unequal variances (Welch's t-test).
inline void tutest(DFMatrixKokkos<double>& data1, int n1,
                   DFMatrixKokkos<double>& data2, int n2,
                   double& t, double& prob)
{
    double ave1, var1, ave2, var2;
    avevar(data1, n1, ave1, var1);
    avevar(data2, n2, ave2, var2);

    t = (ave1 - ave2) / sqrt(var1 / n1 + var2 / n2);
    double df = (var1 / n1 + var2 / n2) * (var1 / n1 + var2 / n2) /
                ((var1 / n1) * (var1 / n1) / (n1 - 1) +
                 (var2 / n2) * (var2 / n2) / (n2 - 1));
    prob = betai(0.5 * df, 0.5, df / (df + t * t));
}
