// Covariance matrix from SVD (Numerical Recipes SVDVAR).
//
// Given the SVD matrices v(1..ma, 1..ma) and singular values w(1..ma),
// computes the covariance matrix cvm(1..ma, 1..ma).
// Singular values that were zeroed (w(k)==0) are excluded from the sum.

#pragma once
#include <matar.h>

using namespace mtr;

inline void svdvar(DFMatrixKokkos<double>& v, int ma,
                   DFMatrixKokkos<double>& w,
                   DFMatrixKokkos<double>& cvm)
{
    double wti[ma];

    for (int i = 1; i <= ma; i++) {
        wti[i - 1] = 0.0;
        if (w.host(i) != 0.0)
            wti[i - 1] = 1.0 / (w.host(i) * w.host(i));
    }

    for (int i = 1; i <= ma; i++) {
        for (int j = 1; j <= i; j++) {
            double sum = 0.0;
            for (int k = 1; k <= ma; k++)
                sum += v.host(i, k) * v.host(j, k) * wti[k - 1];
            cvm.host(i, j) = sum;
            cvm.host(j, i) = sum;
        }
    }
}
