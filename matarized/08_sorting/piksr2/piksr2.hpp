#pragma once
#include <matar.h>

using namespace mtr;

// Insertion sort of two arrays. Sorts arr(1..n) into ascending order while
// making the same rearrangement of brr(1..n).
inline void piksr2(int n, DFMatrixKokkos<double>& arr, DFMatrixKokkos<double>& brr)
{
    for (int j = 2; j <= n; j++) {
        double a = arr.host(j);
        double b = brr.host(j);
        int i;
        for (i = j - 1; i >= 1; i--) {
            if (arr.host(i) <= a) break;
            arr.host(i + 1) = arr.host(i);
            brr.host(i + 1) = brr.host(i);
        }
        arr.host(i + 1) = a;
        brr.host(i + 1) = b;
    }
}
