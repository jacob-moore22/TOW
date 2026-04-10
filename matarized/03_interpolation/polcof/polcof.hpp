// Polynomial coefficients via synthetic division (Numerical Recipes POLCOF).
// Given arrays xa[1..n] and ya[1..n], returns polynomial coefficients
// cof[1..n] such that y = cof(1) + cof(2)*x + ... + cof(n)*x^(n-1).
// Uses polint for intermediate extrapolation to x=0.

#pragma once
#include <cmath>
#include <matar.h>
#include "polint.hpp"

using namespace mtr;

inline void polcof(DFMatrixKokkos<double>& xa, DFMatrixKokkos<double>& ya,
                   int n, DFMatrixKokkos<double>& cof)
{
    constexpr int NMAX = 15;
    DFMatrixKokkos<double> x_loc(NMAX);
    DFMatrixKokkos<double> y_loc(NMAX);

    for (int j = 1; j <= n; j++) {
        x_loc.host(j) = xa.host(j);
        y_loc.host(j) = ya.host(j);
    }

    for (int j = 1; j <= n; j++) {
        double dy;
        polint(x_loc, y_loc, n + 1 - j, 0.0, cof.host(j), dy);
        double xmin = 1.0e38;
        int k = 0;
        for (int i = 1; i <= n + 1 - j; i++) {
            if (std::fabs(x_loc.host(i)) < xmin) {
                xmin = std::fabs(x_loc.host(i));
                k = i;
            }
            if (x_loc.host(i) != 0.0) {
                y_loc.host(i) = (y_loc.host(i) - cof.host(j)) / x_loc.host(i);
            }
        }
        if (k < n + 1 - j) {
            for (int i = k + 1; i <= n + 1 - j; i++) {
                y_loc.host(i - 1) = y_loc.host(i);
                x_loc.host(i - 1) = x_loc.host(i);
            }
        }
    }
}
