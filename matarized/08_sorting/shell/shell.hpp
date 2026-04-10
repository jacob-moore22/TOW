#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

// Shell sort. Sorts arr(1..n) into ascending order using diminishing-increment
// insertion sort (Shell's method).
inline void shell(int n, DFMatrixKokkos<double>& arr)
{
    constexpr double aln2i = 1.4426950;
    constexpr double tiny  = 1.0e-5;

    int lognb2 = static_cast<int>(std::log(static_cast<double>(n)) * aln2i + tiny);
    int m = n;

    for (int nn = 1; nn <= lognb2; nn++) {
        m /= 2;
        int k = n - m;
        for (int j = 1; j <= k; j++) {
            int i = j;
            while (true) {
                int l = i + m;
                if (arr.host(l) < arr.host(i)) {
                    double t  = arr.host(i);
                    arr.host(i) = arr.host(l);
                    arr.host(l) = t;
                    i -= m;
                    if (i < 1) break;
                } else {
                    break;
                }
            }
        }
    }
}
