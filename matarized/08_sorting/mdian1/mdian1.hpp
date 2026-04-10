#pragma once
#include <matar.h>
#include "../sort/sort.hpp"

using namespace mtr;

// Median via sort. Returns the median of x(1..n). The input array x is sorted
// as a side-effect (uses heapsort).
inline void mdian1(DFMatrixKokkos<double>& x, int n, double& xmed)
{
    sort(n, x);
    int n2 = n / 2;
    if (2 * n2 == n) {
        xmed = 0.5 * (x.host(n2) + x.host(n2 + 1));
    } else {
        xmed = x.host(n2 + 1);
    }
}
