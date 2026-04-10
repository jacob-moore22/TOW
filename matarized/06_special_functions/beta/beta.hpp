#pragma once
#include <cmath>
#include <matar.h>
#include "gammln.hpp"

// Beta function B(z,w) = Gamma(z)*Gamma(w) / Gamma(z+w)
// Computed in log space to avoid overflow.
KOKKOS_INLINE_FUNCTION
double beta(double z, double w)
{
    return exp(gammln(z) + gammln(w) - gammln(z + w));
}
