#pragma once
#include <matar.h>

// Extended trapezoidal rule.  Call sequentially with n = 1, 2, 3, ...
// On the first call (n == 1) it computes the crudest estimate; each
// subsequent call refines the previous estimate stored in s.
// The caller must keep s and it alive between calls.
template<typename Func>
KOKKOS_INLINE_FUNCTION
void trapzd(Func func, double a, double b, double& s, int n, int& it)
{
    if (n == 1) {
        s  = 0.5 * (b - a) * (func(a) + func(b));
        it = 1;
    } else {
        double tnm = it;
        double del = (b - a) / tnm;
        double x   = a + 0.5 * del;
        double sum = 0.0;
        for (int j = 0; j < it; j++) {
            sum += func(x);
            x   += del;
        }
        s   = 0.5 * (s + (b - a) * sum / tnm);
        it *= 2;
    }
}
