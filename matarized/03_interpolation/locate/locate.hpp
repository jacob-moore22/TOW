// Binary search in an ordered table (Numerical Recipes LOCATE).
// Given array xx[1..n] and a value x, returns j such that x is
// between xx(j) and xx(j+1). j=0 or j=n indicates out of range.

#pragma once
#include <matar.h>

using namespace mtr;

inline void locate(DFMatrixKokkos<double>& xx, int n, double x, int& j)
{
    int jl = 0;
    int ju = n + 1;
    while (ju - jl > 1) {
        int jm = (ju + jl) / 2;
        if ((xx.host(n) > xx.host(1)) == (x > xx.host(jm))) {
            jl = jm;
        } else {
            ju = jm;
        }
    }
    j = jl;
}
