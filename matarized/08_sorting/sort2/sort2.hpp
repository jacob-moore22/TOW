#pragma once
#include <matar.h>

using namespace mtr;

// Heapsort of two arrays. Sorts ra(1..n) into ascending order using heapsort
// while making the same rearrangement of rb(1..n).
inline void sort2(int n, DFMatrixKokkos<double>& ra, DFMatrixKokkos<double>& rb)
{
    if (n < 2) return;

    int l  = n / 2 + 1;
    int ir = n;

    while (true) {
        double rra, rrb;
        if (l > 1) {
            l--;
            rra = ra.host(l);
            rrb = rb.host(l);
        } else {
            rra = ra.host(ir);
            rrb = rb.host(ir);
            ra.host(ir) = ra.host(1);
            rb.host(ir) = rb.host(1);
            ir--;
            if (ir == 1) {
                ra.host(1) = rra;
                rb.host(1) = rrb;
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
                rb.host(i) = rb.host(j);
                i = j;
                j = j + j;
            } else {
                j = ir + 1;
            }
        }
        ra.host(i) = rra;
        rb.host(i) = rrb;
    }
}
