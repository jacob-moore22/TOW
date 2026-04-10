// Bicubic interpolation coefficients (Numerical Recipes BCUCOF).
// Given function values y, gradients y1 (d/dx1), y2 (d/dx2), and
// cross-derivatives y12 at the four corners of a grid cell, plus
// the cell dimensions d1 and d2, returns the 4x4 coefficient matrix c.
//
// All four-element arrays are indexed 1..4 (corners: ll, lr, ur, ul).
// c is a 4x4 DFMatrixKokkos (1-based).

#pragma once
#include <matar.h>

using namespace mtr;

inline void bcucof(DFMatrixKokkos<double>& y, DFMatrixKokkos<double>& y1,
                   DFMatrixKokkos<double>& y2, DFMatrixKokkos<double>& y12,
                   double d1, double d2,
                   DFMatrixKokkos<double>& c)
{
    // Weight matrix stored in Fortran column-major order (matches DATA statement)
    static constexpr double wt_flat[256] = {
         1, 0,-3, 2, 0, 0, 0, 0,-3, 0, 9,-6, 2, 0,-6, 4,
         0, 0, 0, 0, 0, 0, 0, 0, 3, 0,-9, 6,-2, 0, 6,-4,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-6, 0, 0,-6, 4,
         0, 0, 3,-2, 0, 0, 0, 0, 0, 0,-9, 6, 0, 0, 6,-4,
         0, 0, 0, 0, 1, 0,-3, 2,-2, 0, 6,-4, 1, 0,-3, 2,
         0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 3,-2, 1, 0,-3, 2,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 2, 0, 0, 3,-2,
         0, 0, 0, 0, 0, 0, 3,-2, 0, 0,-6, 4, 0, 0, 3,-2,
         0, 1,-2, 1, 0, 0, 0, 0, 0,-3, 6,-3, 0, 2,-4, 2,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 3,-6, 3, 0,-2, 4,-2,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 2,-2,
         0, 0,-1, 1, 0, 0, 0, 0, 0, 0, 3,-3, 0, 0,-2, 2,
         0, 0, 0, 0, 0, 1,-2, 1, 0,-2, 4,-2, 0, 1,-2, 1,
         0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 2,-1, 0, 1,-2, 1,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1, 0, 0,-1, 1,
         0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 2,-2, 0, 0,-1, 1
    };

    double d1d2 = d1 * d2;
    double x_vec[17]; // 1-based
    for (int i = 1; i <= 4; i++) {
        x_vec[i]      = y.host(i);
        x_vec[i + 4]  = y1.host(i) * d1;
        x_vec[i + 8]  = y2.host(i) * d2;
        x_vec[i + 12] = y12.host(i) * d1d2;
    }

    double cl[17]; // 1-based
    for (int i = 1; i <= 16; i++) {
        double xx = 0.0;
        for (int k = 1; k <= 16; k++) {
            // Fortran WT(I,K) = wt_flat[(K-1)*16 + (I-1)]
            xx += wt_flat[(k - 1) * 16 + (i - 1)] * x_vec[k];
        }
        cl[i] = xx;
    }

    int l = 0;
    for (int i = 1; i <= 4; i++) {
        for (int j = 1; j <= 4; j++) {
            l++;
            c.host(i, j) = cl[l];
        }
    }
}
