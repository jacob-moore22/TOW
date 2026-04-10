#pragma once
#include <cmath>
#include <matar.h>
#include "gammp.hpp"
#include "gammq.hpp"

// Complementary error function erfc(x) = 1 - erf(x).
// Named nr_erfc to avoid collision with std::erfc.
KOKKOS_INLINE_FUNCTION
double nr_erfc(double x)
{
    if (x < 0.0) {
        return 1.0 + gammp(0.5, x * x);
    } else {
        return gammq(0.5, x * x);
    }
}
