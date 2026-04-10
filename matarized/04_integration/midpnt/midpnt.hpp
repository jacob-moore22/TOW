#pragma once
#include <matar.h>

// Open-formula midpoint rule.  Call sequentially with n = 1, 2, 3, ...
// On n == 1 it returns the crudest estimate; subsequent calls triple the
// number of interior points.  Caller keeps s and it alive between calls.
template<typename Func>
KOKKOS_INLINE_FUNCTION
void midpnt(Func func, double a, double b, double& s, int n, int& it)
{
    if (n == 1) {
        s  = (b - a) * func(0.5 * (a + b));
        it = 1;
    } else {
        double tnm  = it;
        double del  = (b - a) / (3.0 * tnm);
        double ddel = del + del;
        double x    = a + 0.5 * del;
        double sum  = 0.0;
        for (int j = 0; j < it; j++) {
            sum += func(x);
            x   += ddel;
            sum += func(x);
            x   += del;
        }
        s   = (s + (b - a) * sum / tnm) / 3.0;
        it *= 3;
    }
}
