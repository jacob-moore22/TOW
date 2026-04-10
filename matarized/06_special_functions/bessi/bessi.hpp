#pragma once
#include <cmath>
#include <matar.h>
#include "bessi0.hpp"
#include "bessi1.hpp"

// Modified Bessel function In(x) for integer n >= 2 via downward recurrence.
// Numerical Recipes (Press et al.), adapted for MATAR/Kokkos.
KOKKOS_INLINE_FUNCTION
double bessi(int n, double x)
{
    constexpr int    IACC  = 40;
    constexpr double BIGNO = 1.0e10;
    constexpr double BIGNI = 1.0e-10;

    if (n < 2) return 0.0;
    if (x == 0.0) return 0.0;

    double tox = 2.0 / x;
    double bip = 0.0;
    double bi  = 1.0;
    double ans = 0.0;

    int m = 2 * (n + static_cast<int>(sqrt(static_cast<double>(IACC * n))));

    for (int j = m; j >= 1; j--) {
        double bim = bip + static_cast<double>(j) * tox * bi;
        bip = bi;
        bi  = bim;
        if (fabs(bi) > BIGNO) {
            ans *= BIGNI;
            bi  *= BIGNI;
            bip *= BIGNI;
        }
        if (j == n) ans = bip;
    }

    return ans * bessi0(x) / bi;
}
