#pragma once
#include <matar.h>

using namespace mtr;

// Reduction step for solvde (Numerical Recipes RED).
// Combines information from the previous mesh point (stored in c)
// with the current set of equations (in s).
inline void red(int iz1, int iz2, int jz1, int jz2,
                int jm1, int jm2, int jmf,
                int ic1, int jc1, int jcf, int kc,
                DFMatrixKokkos<double>& c, int /*nci*/, int /*ncj*/, int /*nck*/,
                DFMatrixKokkos<double>& s, int /*nsi*/, int /*nsj*/)
{
    int loff = jc1 - jm1;
    int ic   = ic1;
    for (int j = jz1; j <= jz2; j++) {
        for (int l = jm1; l <= jm2; l++) {
            double vx = c.host(ic, l + loff, kc);
            for (int i = iz1; i <= iz2; i++)
                s.host(i, l) -= s.host(i, j) * vx;
        }
        double vx = c.host(ic, jcf, kc);
        for (int i = iz1; i <= iz2; i++)
            s.host(i, jmf) -= s.host(i, j) * vx;
        ic++;
    }
}
