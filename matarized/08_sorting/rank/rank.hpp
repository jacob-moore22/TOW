#pragma once
#include <matar.h>

using namespace mtr;

// Given an index array indx(1..n) as output by indexx, produces a rank table
// irank(1..n) such that irank(j) is the rank of the j-th element of the
// original array.
inline void rank(int n, DFMatrixKokkos<int>& indx, DFMatrixKokkos<int>& irank)
{
    for (int j = 1; j <= n; j++) {
        irank.host(indx.host(j)) = j;
    }
}
