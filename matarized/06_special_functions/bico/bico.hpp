#pragma once
#include <cmath>
#include <matar.h>
#include "factln.hpp"

// Binomial coefficient C(n,k) = n! / (k! * (n-k)!)
// Computed in log space to avoid overflow, then rounded.
KOKKOS_INLINE_FUNCTION
double bico(int n, int k)
{
    return round(exp(factln(n) - factln(k) - factln(n - k)));
}
