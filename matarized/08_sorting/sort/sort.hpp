#pragma once
#include <matar.h>

using namespace mtr;

// Heapsort. Sorts ra(1..n) into ascending order using the heapsort algorithm.
inline void sort(int n, DFMatrixKokkos<double>& ra)
{
    if (n < 2) return;

    int l  = n / 2 + 1;
    int ir = n;

    while (true) {
        double rra;
        if (l > 1) {
            l--;
            rra = ra.host(l);
        } else {
            rra = ra.host(ir);
            ra.host(ir) = ra.host(1);
            ir--;
            if (ir == 1) {
                ra.host(1) = rra;
                return;
            }
        }
        int i = l;
        int j = l + l;
        while (j <= ir) {
            if (j < ir) {
                if (ra.host(j) < ra.host(j + 1)) j++;
            }
            if (rra < ra.host(j)) {
                ra.host(i) = ra.host(j);
                i = j;
                j = j + j;
            } else {
                j = ir + 1;
            }
        }
        ra.host(i) = rra;
    }
}
