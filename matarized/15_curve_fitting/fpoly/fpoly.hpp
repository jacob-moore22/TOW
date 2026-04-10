// Polynomial basis functions (Numerical Recipes FPOLY).
//
// Evaluates p(1..np) = {1, x, x^2, ..., x^(np-1)} for use with lfit/svdfit.

#pragma once
#include <matar.h>

using namespace mtr;

inline void fpoly(double x, DFMatrixKokkos<double>& p, int np)
{
    p.host(1) = 1.0;
    for (int j = 2; j <= np; j++)
        p.host(j) = p.host(j - 1) * x;
}
