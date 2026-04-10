#pragma once
#include <matar.h>
#include "../indexx/indexx.hpp"

using namespace mtr;

// Index-based sort of three arrays. Sorts ra(1..n) into ascending order using
// indexx, making the same rearrangement of rb(1..n) and rc(1..n). Uses wksp
// as a temporary workspace and iwksp for the index array.
inline void sort3(int n,
                  DFMatrixKokkos<double>& ra,
                  DFMatrixKokkos<double>& rb,
                  DFMatrixKokkos<double>& rc,
                  DFMatrixKokkos<double>& wksp,
                  DFMatrixKokkos<int>&    iwksp)
{
    indexx(n, ra, iwksp);

    for (int j = 1; j <= n; j++) {
        wksp.host(j) = ra.host(j);
    }
    for (int j = 1; j <= n; j++) {
        ra.host(j) = wksp.host(iwksp.host(j));
    }

    for (int j = 1; j <= n; j++) {
        wksp.host(j) = rb.host(j);
    }
    for (int j = 1; j <= n; j++) {
        rb.host(j) = wksp.host(iwksp.host(j));
    }

    for (int j = 1; j <= n; j++) {
        wksp.host(j) = rc.host(j);
    }
    for (int j = 1; j <= n; j++) {
        rc.host(j) = wksp.host(iwksp.host(j));
    }
}
