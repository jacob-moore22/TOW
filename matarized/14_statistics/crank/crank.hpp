#pragma once
#include <matar.h>

using namespace mtr;

// Replace sorted array w(1..n) by its ranks, handling ties.
// Returns s = sum of (t^3 - t) for groups of t tied values.
inline void crank(int n, DFMatrixKokkos<double>& w, double& s)
{
    s = 0.0;
    int j = 1;
    while (j < n) {
        if (w.host(j + 1) != w.host(j)) {
            w.host(j) = static_cast<double>(j);
            j++;
        } else {
            int jt;
            for (jt = j + 1; jt <= n; jt++) {
                if (w.host(jt) != w.host(j)) break;
            }
            double rank = 0.5 * (j + jt - 1);
            for (int ji = j; ji <= jt - 1; ji++)
                w.host(ji) = rank;
            double t = jt - j;
            s += t * t * t - t;
            j = jt;
        }
    }
    if (j == n) w.host(n) = static_cast<double>(n);
}
