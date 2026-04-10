#pragma once
#include <matar.h>

using namespace mtr;

// Back-substitution for the block-banded matrix from solvde
// (Numerical Recipes BKSUB).
// c is a 3D array c(nci, ncj, nck) with 1-based Fortran indexing.
inline void bksub(int ne, int nb, int jf, int k1, int k2,
                  DFMatrixKokkos<double>& c, int nci, int ncj, int nck)
{
    int nbf = ne - nb;

    for (int k = k2; k >= k1; k--) {
        int kp = k + 1;
        for (int j = 1; j <= nbf; j++) {
            double xx = c.host(j, jf, kp);
            for (int i = 1; i <= ne; i++)
                c.host(i, jf, k) = c.host(i, jf, k) - c.host(i, j, k) * xx;
        }
    }

    for (int k = k1; k <= k2; k++) {
        int kp = k + 1;
        for (int i = 1; i <= nb; i++)
            c.host(i, 1, k) = c.host(i + nbf, jf, k);
        for (int i = 1; i <= nbf; i++)
            c.host(i + nb, 1, k) = c.host(i, jf, kp);
    }
}
