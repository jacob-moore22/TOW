#pragma once
#include <cmath>
#include <matar.h>

// Associated Legendre polynomial P_l^m(x) via the standard recurrence.
// Requires 0 <= m <= l and |x| <= 1.
KOKKOS_INLINE_FUNCTION
double plgndr(int l, int m, double x)
{
    double pmm = 1.0;
    if (m > 0) {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact  = 1.0;
        for (int i = 1; i <= m; i++) {
            pmm  = -pmm * fact * somx2;
            fact += 2.0;
        }
    }

    if (l == m) {
        return pmm;
    }

    double pmmp1 = x * (2 * m + 1) * pmm;
    if (l == m + 1) {
        return pmmp1;
    }

    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ll++) {
        pll   = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
        pmm   = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}
