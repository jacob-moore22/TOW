#pragma once
#include <cmath>
#include <matar.h>

// Log of the gamma function using the Lanczos approximation.
// KOKKOS_INLINE_FUNCTION allows calling from within FOR_ALL / DO_ALL kernels.
KOKKOS_INLINE_FUNCTION
double gammln(double xx)
{
    constexpr double cof[6] = {
         76.18009173,
        -86.50532033,
         24.01409822,
         -1.231739516,
          0.120858003e-2,
         -0.536382e-5
    };
    constexpr double stp = 2.50662827465;

    double x   = xx - 1.0;
    double tmp = x + 5.5;
    tmp = (x + 0.5) * log(tmp) - tmp;

    double ser = 1.0;
    for (int j = 0; j < 6; j++) {
        x   += 1.0;
        ser += cof[j] / x;
    }
    return tmp + log(stp * ser);
}
