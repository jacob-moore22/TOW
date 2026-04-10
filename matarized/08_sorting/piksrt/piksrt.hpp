#pragma once
#include <matar.h>

using namespace mtr;

// Insertion sort. Sorts arr(1..n) into ascending order.
inline void piksrt(int n, DFMatrixKokkos<double>& arr)
{
    for (int j = 2; j <= n; j++) {
        double a = arr.host(j);
        int i;
        for (i = j - 1; i >= 1; i--) {
            if (arr.host(i) <= a) break;
            arr.host(i + 1) = arr.host(i);
        }
        arr.host(i + 1) = a;
    }
}
