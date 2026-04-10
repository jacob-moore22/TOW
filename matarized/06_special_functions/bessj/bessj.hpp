#pragma once
#include <cmath>
#include <matar.h>
#include "bessj0.hpp"
#include "bessj1.hpp"

// Bessel function Jn(x) for integer n >= 2.
// Uses upward recurrence when x > n, Miller's downward recurrence otherwise.
// Numerical Recipes (Press et al.)
KOKKOS_INLINE_FUNCTION
double bessj(int n, double x)
{
    constexpr int    IACC  = 40;
    constexpr double BIGNO = 1.0e10;
    constexpr double BIGNI = 1.0e-10;

    if (n < 2) return 0.0;

    double tox = 2.0 / x;

    if (x > static_cast<double>(n)) {
        double bjm = bessj0(x);
        double bj  = bessj1(x);
        for (int j = 1; j <= n - 1; j++) {
            double bjp = j * tox * bj - bjm;
            bjm = bj;
            bj  = bjp;
        }
        return bj;
    } else {
        int m = 2 * ((n + static_cast<int>(sqrt(static_cast<double>(IACC * n)))) / 2);
        double result = 0.0;
        int    jsum   = 0;
        double sum    = 0.0;
        double bjp    = 0.0;
        double bj     = 1.0;

        for (int j = m; j >= 1; j--) {
            double bjm = j * tox * bj - bjp;
            bjp = bj;
            bj  = bjm;
            if (fabs(bj) > BIGNO) {
                bj     *= BIGNI;
                bjp    *= BIGNI;
                result *= BIGNI;
                sum    *= BIGNI;
            }
            if (jsum != 0) sum += bj;
            jsum = 1 - jsum;
            if (j == n) result = bjp;
        }
        sum = 2.0 * sum - bj;
        return result / sum;
    }
}
