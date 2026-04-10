#pragma once
#include <matar.h>

// Open midpoint rule for semi-infinite intervals.
// Evaluates the integral of funk over [aa, bb] where one or both limits
// may be very large (or infinite).  Internally maps [aa, bb] -> [1/bb, 1/aa]
// via the substitution x = 1/t.
// Call sequentially with n = 1, 2, 3, ... (same protocol as midpnt).
template<typename Func>
KOKKOS_INLINE_FUNCTION
void midinf(Func funk, double aa, double bb, double& s, int n, int& it)
{
    double b = 1.0 / aa;
    double a = 1.0 / bb;

    if (n == 1) {
        double xmid = 0.5 * (a + b);
        s  = (b - a) * funk(1.0 / xmid) / (xmid * xmid);
        it = 1;
    } else {
        double tnm  = it;
        double del  = (b - a) / (3.0 * tnm);
        double ddel = del + del;
        double x    = a + 0.5 * del;
        double sum  = 0.0;
        for (int j = 0; j < it; j++) {
            sum += funk(1.0 / x) / (x * x);
            x   += ddel;
            sum += funk(1.0 / x) / (x * x);
            x   += del;
        }
        s   = (s + (b - a) * sum / tnm) / 3.0;
        it *= 3;
    }
}
