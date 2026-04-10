#pragma once
#include <matar.h>

using namespace mtr;

// Equivalence classes from pairs. Given m pairs of equivalent elements in
// lista(1..m) and listb(1..m), determines the equivalence classes for n
// elements, returning class labels in nf(1..n).
inline void eclass(DFMatrixKokkos<int>& nf, int n,
                   DFMatrixKokkos<int>& lista, DFMatrixKokkos<int>& listb, int m)
{
    for (int k = 1; k <= n; k++) {
        nf.host(k) = k;
    }

    for (int l = 1; l <= m; l++) {
        int j = lista.host(l);
        while (nf.host(j) != j) {
            j = nf.host(j);
        }
        int k = listb.host(l);
        while (nf.host(k) != k) {
            k = nf.host(k);
        }
        if (j != k) nf.host(j) = k;
    }

    for (int j = 1; j <= n; j++) {
        while (nf.host(j) != nf.host(nf.host(j))) {
            nf.host(j) = nf.host(nf.host(j));
        }
    }
}
