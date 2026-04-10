#pragma once
#include <cmath>
#include <matar.h>
#include "erfcc.hpp"

using namespace mtr;

// Kendall's tau for a contingency table tab(1..ni, 1..nj).
// Returns tau, z, and significance prob.
inline void kendl2(DFMatrixKokkos<double>& tab, int ni, int nj,
                   double& tau, double& z, double& prob)
{
    double en1 = 0.0, en2 = 0.0, s = 0.0;
    int nn = ni * nj;

    double points = tab.host(ni, nj);
    for (int k = 0; k <= nn - 2; k++) {
        int ki = k / nj;
        int kj = k - nj * ki;
        points += tab.host(ki + 1, kj + 1);
        for (int l = k + 1; l <= nn - 1; l++) {
            int li = l / nj;
            int lj = l - nj * li;
            int m1 = li - ki;
            int m2 = lj - kj;
            int mm = m1 * m2;
            double pairs = tab.host(ki + 1, kj + 1) * tab.host(li + 1, lj + 1);
            if (mm != 0) {
                en1 += pairs;
                en2 += pairs;
                if (mm > 0) s += pairs; else s -= pairs;
            } else {
                if (m1 != 0) en1 += pairs;
                if (m2 != 0) en2 += pairs;
            }
        }
    }

    tau = s / sqrt(en1 * en2);
    double var = (4.0 * points + 10.0) / (9.0 * points * (points - 1.0));
    z = tau / sqrt(var);
    prob = erfcc(fabs(z) / 1.4142136);
}
