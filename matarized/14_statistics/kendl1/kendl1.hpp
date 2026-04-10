#pragma once
#include <cmath>
#include <matar.h>
#include "erfcc.hpp"

using namespace mtr;

// Kendall's tau rank correlation for two vectors data1(1..n) and data2(1..n).
// Returns tau, its normal deviate z, and significance prob.
inline void kendl1(DFMatrixKokkos<double>& data1, DFMatrixKokkos<double>& data2,
                   int n, double& tau, double& z, double& prob)
{
    long n1 = 0, n2 = 0;
    long is = 0;

    for (int j = 1; j <= n - 1; j++) {
        for (int k = j + 1; k <= n; k++) {
            double a1 = data1.host(j) - data1.host(k);
            double a2 = data2.host(j) - data2.host(k);
            double aa = a1 * a2;
            if (aa != 0.0) {
                n1++;
                n2++;
                if (aa > 0.0) is++; else is--;
            } else {
                if (a1 != 0.0) n1++;
                if (a2 != 0.0) n2++;
            }
        }
    }

    tau = static_cast<double>(is) / sqrt(static_cast<double>(n1) * static_cast<double>(n2));
    double var = (4.0 * n + 10.0) / (9.0 * n * (n - 1.0));
    z = tau / sqrt(var);
    prob = erfcc(fabs(z) / 1.4142136);
}
