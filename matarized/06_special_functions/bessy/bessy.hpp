#pragma once
#include <cmath>
#include <matar.h>
#include "bessy0.hpp"
#include "bessy1.hpp"

// Bessel function Yn(x) for integer n >= 2, x > 0.
// Uses upward recurrence from Y0 and Y1.
// Numerical Recipes (Press et al.)
KOKKOS_INLINE_FUNCTION
double bessy(int n, double x)
{
    if (n < 2) return 0.0;

    double tox = 2.0 / x;
    double by  = bessy1(x);
    double bym = bessy0(x);

    for (int j = 1; j <= n - 1; j++) {
        double byp = j * tox * by - bym;
        bym = by;
        by  = byp;
    }
    return by;
}
