#pragma once
#include <cmath>
#include <matar.h>
#include "gammp.hpp"

// Error function erf(x). Named nr_erf to avoid collision with std::erf.
KOKKOS_INLINE_FUNCTION
double nr_erf(double x)
{
    if (x < 0.0) {
        return -gammp(0.5, x * x);
    } else {
        return  gammp(0.5, x * x);
    }
}
