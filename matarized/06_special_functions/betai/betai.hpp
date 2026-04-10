#pragma once
#include <cmath>
#include <matar.h>
#include "gammln.hpp"
#include "betacf.hpp"

// Incomplete beta function I_x(a,b).
KOKKOS_INLINE_FUNCTION
double betai(double a, double b, double x)
{
    double bt;
    if (x == 0.0 || x == 1.0) {
        bt = 0.0;
    } else {
        bt = exp(gammln(a + b) - gammln(a) - gammln(b)
             + a * log(x) + b * log(1.0 - x));
    }

    if (x < (a + 1.0) / (a + b + 2.0)) {
        return bt * betacf(a, b, x) / a;
    } else {
        return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
    }
}
