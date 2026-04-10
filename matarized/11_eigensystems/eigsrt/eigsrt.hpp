// Sort eigenvalues and eigenvectors in descending order (Numerical Recipes EIGSRT).
// Given eigenvalues d[1..n] and eigenvectors v[1..n][1..n] (columns),
// rearranges them so that d[1] >= d[2] >= ... >= d[n].

#pragma once
#include <matar.h>

using namespace mtr;

inline void eigsrt(DFMatrixKokkos<double>& d, DFMatrixKokkos<double>& v, int n)
{
    for (int i = 1; i <= n - 1; i++) {
        int k = i;
        double p = d.host(i);
        for (int j = i + 1; j <= n; j++) {
            if (d.host(j) >= p) {
                k = j;
                p = d.host(j);
            }
        }
        if (k != i) {
            d.host(k) = d.host(i);
            d.host(i) = p;
            for (int j = 1; j <= n; j++) {
                p = v.host(j, i);
                v.host(j, i) = v.host(j, k);
                v.host(j, k) = p;
            }
        }
    }
}
