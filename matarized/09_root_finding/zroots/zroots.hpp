#pragma once
#include <cmath>
#include <complex>
#include <cstdio>
#include <matar.h>
#include "../laguer/laguer.hpp"

using namespace mtr;

// Find all m roots of a polynomial with complex coefficients a[0..m].
// Uses successive deflation via Laguerre's method, then optionally polishes
// each root against the original (undeflated) polynomial. Roots are sorted
// by real part on exit.
inline void zroots(const std::complex<double>* a, int m,
                   std::complex<double>* roots, bool polish)
{
    constexpr double EPS  = 1.0e-6;
    constexpr int    MAXM = 101;

    std::complex<double> ad[MAXM];

    for (int j = 0; j <= m; j++)
        ad[j] = a[j];

    for (int j = m; j >= 1; j--) {
        std::complex<double> x(0.0, 0.0);
        laguer(ad, j, x, EPS, false);

        if (std::fabs(x.imag()) <= 2.0 * EPS * EPS * std::fabs(x.real()))
            x = std::complex<double>(x.real(), 0.0);

        roots[j - 1] = x;

        std::complex<double> b = ad[j];
        for (int jj = j - 1; jj >= 0; jj--) {
            std::complex<double> c = ad[jj];
            ad[jj] = b;
            b = x * b + c;
        }
    }

    if (polish) {
        for (int j = 0; j < m; j++)
            laguer(a, m, roots[j], EPS, true);
    }

    // Insertion sort by real part
    for (int j = 1; j < m; j++) {
        std::complex<double> x = roots[j];
        int i;
        for (i = j - 1; i >= 0; i--) {
            if (roots[i].real() <= x.real()) break;
            roots[i + 1] = roots[i];
        }
        roots[i + 1] = x;
    }
}
