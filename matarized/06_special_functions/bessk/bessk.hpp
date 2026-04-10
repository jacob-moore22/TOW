#pragma once
#include <cmath>
#include <matar.h>
#include "bessk0.hpp"
#include "bessk1.hpp"

// Modified Bessel function Kn(x) for integer n >= 2 via upward recurrence.
// Numerical Recipes (Press et al.), adapted for MATAR/Kokkos.
KOKKOS_INLINE_FUNCTION
double bessk(int n, double x)
{
    if (n < 2) return 0.0;

    double tox = 2.0 / x;
    double bkm = bessk0(x);
    double bk  = bessk1(x);

    for (int j = 1; j <= n - 1; j++) {
        double bkp = bkm + static_cast<double>(j) * tox * bk;
        bkm = bk;
        bk  = bkp;
    }

    return bk;
}
