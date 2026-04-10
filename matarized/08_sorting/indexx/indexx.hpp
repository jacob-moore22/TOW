#pragma once
#include <matar.h>

using namespace mtr;

// Indirect heapsort. Indexes arrin(1..n), outputting the index array indx(1..n)
// such that arrin(indx(j)) is in ascending order for j = 1..n. arrin is not
// modified.
inline void indexx(int n, DFMatrixKokkos<double>& arrin, DFMatrixKokkos<int>& indx)
{
    for (int j = 1; j <= n; j++) {
        indx.host(j) = j;
    }

    if (n < 2) return;

    int l  = n / 2 + 1;
    int ir = n;

    while (true) {
        int indxt;
        double q;
        if (l > 1) {
            l--;
            indxt = indx.host(l);
            q = arrin.host(indxt);
        } else {
            indxt = indx.host(ir);
            q = arrin.host(indxt);
            indx.host(ir) = indx.host(1);
            ir--;
            if (ir == 1) {
                indx.host(1) = indxt;
                return;
            }
        }
        int i = l;
        int j = l + l;
        while (j <= ir) {
            if (j < ir) {
                if (arrin.host(indx.host(j)) < arrin.host(indx.host(j + 1))) j++;
            }
            if (q < arrin.host(indx.host(j))) {
                indx.host(i) = indx.host(j);
                i = j;
                j = j + j;
            } else {
                j = ir + 1;
            }
        }
        indx.host(i) = indxt;
    }
}
