#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

// Helper for the simplex LP solver. Finds the maximum element in row mm of the
// tableau among columns indexed by ll[0..nll-1]. If iabf != 0, compares by
// absolute value. Returns the column index kp and the value bmax.
//
// a:    tableau stored row-major, dimensions (mp rows) x (np cols).
// np:   number of columns in the tableau.
// mm:   row to scan (0-based).
// ll:   array of column indices to consider (0-based), length nll.
// nll:  number of entries in ll.
// iabf: if nonzero, compare absolute values.
// kp:   on return, column index of the maximum.
// bmax: on return, the maximum value found.
inline void simp1(const double* a, int np, int mm, const int* ll, int nll,
                  int iabf, int& kp, double& bmax)
{
    kp   = ll[0];
    bmax = a[mm * np + kp];
    if (nll < 2) return;

    for (int k = 1; k < nll; k++) {
        double test;
        if (iabf == 0)
            test = a[mm * np + ll[k]] - bmax;
        else
            test = std::fabs(a[mm * np + ll[k]]) - std::fabs(bmax);

        if (test > 0.0) {
            bmax = a[mm * np + ll[k]];
            kp   = ll[k];
        }
    }
}
