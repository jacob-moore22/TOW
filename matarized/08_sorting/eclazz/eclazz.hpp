#pragma once
#include <matar.h>

using namespace mtr;

// Equivalence classes from a function. Given an equivalence-testing function
// equiv(i, j) that returns true if elements i and j are equivalent, determines
// the equivalence classes for n elements, returning class labels in nf(1..n).
template <typename EquivFunc>
inline void eclazz(DFMatrixKokkos<int>& nf, int n, EquivFunc equiv)
{
    nf.host(1) = 1;
    for (int jj = 2; jj <= n; jj++) {
        nf.host(jj) = jj;
        for (int kk = 1; kk <= jj - 1; kk++) {
            nf.host(kk) = nf.host(nf.host(kk));
            if (equiv(jj, kk)) {
                nf.host(nf.host(nf.host(kk))) = jj;
            }
        }
    }
    for (int jj = 1; jj <= n; jj++) {
        nf.host(jj) = nf.host(nf.host(jj));
    }
}
