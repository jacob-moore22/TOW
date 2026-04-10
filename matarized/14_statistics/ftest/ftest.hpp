#pragma once
#include <cmath>
#include <matar.h>
#include "avevar.hpp"
#include "betai.hpp"

using namespace mtr;

// F-test for difference in variances between data1(1..n1) and data2(1..n2).
// Returns F statistic and probability prob.
inline void ftest(DFMatrixKokkos<double>& data1, int n1,
                  DFMatrixKokkos<double>& data2, int n2,
                  double& f, double& prob)
{
    double ave1, var1, ave2, var2;
    avevar(data1, n1, ave1, var1);
    avevar(data2, n2, ave2, var2);

    double df1, df2;
    if (var1 > var2) {
        f   = var1 / var2;
        df1 = n1 - 1;
        df2 = n2 - 1;
    } else {
        f   = var2 / var1;
        df1 = n2 - 1;
        df2 = n1 - 1;
    }

    prob = betai(0.5 * df2, 0.5 * df1, df2 / (df2 + df1 * f))
         + (1.0 - betai(0.5 * df1, 0.5 * df2, df1 / (df1 + df2 / f)));
}
