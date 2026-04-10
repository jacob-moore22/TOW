#pragma once
#include <matar.h>

using namespace mtr;

// Helper for the simplex LP solver. Performs a pivot operation on the tableau,
// eliminating column kp using row ip as the pivot row.
//
// a:    tableau, row-major, dimensions at least (i1+2) rows x (k1+2) cols.
// np:   column stride of the tableau.
// i1:   index of the last row to update (0-based; objective row is 0).
// k1:   index of the last variable column (0-based).
// ip:   pivot row (0-based).
// kp:   pivot column (0-based).
inline void simp3(double* a, int np, int i1, int k1, int ip, int kp)
{
    double piv = 1.0 / a[ip * np + kp];

    if (i1 >= 0) {
        for (int ii = 0; ii <= i1; ii++) {
            if (ii != ip) {
                a[ii * np + kp] *= piv;
                for (int kk = 0; kk <= k1; kk++) {
                    if (kk != kp)
                        a[ii * np + kk] -= a[ip * np + kk] * a[ii * np + kp];
                }
            }
        }
    }

    for (int kk = 0; kk <= k1; kk++) {
        if (kk != kp) a[ip * np + kk] *= -piv;
    }
    a[ip * np + kp] = piv;
}
