#pragma once
#include <cmath>
#include <matar.h>
#include "gammln.hpp"

// Log factorial ln(n!) = ln(Gamma(n+1)).
// The original Fortran caches results for n <= 99;
// device kernels cannot use static tables, so we compute directly.
KOKKOS_INLINE_FUNCTION
double factln(int n)
{
    if (n < 0) {
        return -1.0;
    }
    if (n <= 1) {
        return 0.0;
    }
    return gammln(n + 1.0);
}
