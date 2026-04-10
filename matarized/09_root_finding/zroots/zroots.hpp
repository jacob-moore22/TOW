// Find all m roots of a polynomial with complex coefficients.
// Uses successive deflation via Laguerre's method, then optionally polishes
// each root against the original (undeflated) polynomial. Roots are sorted
// by real part on exit.
//
// a:      interleaved real/imag coefficients, DFMatrixKokkos<double>(2*(m+1)).
//         Coefficient j (0-based): real at a(2*j+1), imag at a(2*j+2).
// m:      polynomial degree
// roots:  interleaved real/imag output, DFMatrixKokkos<double>(2*m).
//         Root j (0-based): real at roots(2*j+1), imag at roots(2*j+2).
// polish: if true, polish each root against the original polynomial.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>
#include "../laguer/laguer.hpp"

using namespace mtr;

inline void zroots(DFMatrixKokkos<double>& a, int m,
                   DFMatrixKokkos<double>& roots, bool polish)
{
    constexpr double EPS  = 1.0e-6;
    constexpr int    MAXM = 101;

    // Working copy of coefficients for deflation (interleaved re/im)
    double ad_re[MAXM], ad_im[MAXM];
    for (int j = 0; j <= m; j++) {
        ad_re[j] = a.host(2 * j + 1);
        ad_im[j] = a.host(2 * j + 2);
    }

    for (int j = m; j >= 1; j--) {
        // Build interleaved DFMatrixKokkos for this deflated polynomial
        DFMatrixKokkos<double> ad_matar(2 * (j + 1));
        for (int k = 0; k <= j; k++) {
            ad_matar.host(2 * k + 1) = ad_re[k];
            ad_matar.host(2 * k + 2) = ad_im[k];
        }

        double xr = 0.0, xi = 0.0;
        laguer(ad_matar, j, xr, xi, EPS, false);

        if (fabs(xi) <= 2.0 * EPS * EPS * fabs(xr))
            xi = 0.0;

        roots.host(2 * (j - 1) + 1) = xr;
        roots.host(2 * (j - 1) + 2) = xi;

        // Deflate: synthetic division by (z - root)
        PortableComplex x(xr, xi);
        PortableComplex b(ad_re[j], ad_im[j]);
        for (int jj = j - 1; jj >= 0; jj--) {
            PortableComplex c(ad_re[jj], ad_im[jj]);
            ad_re[jj] = b.re;
            ad_im[jj] = b.im;
            b = x * b + c;
        }
    }

    if (polish) {
        for (int j = 0; j < m; j++) {
            double xr = roots.host(2 * j + 1);
            double xi = roots.host(2 * j + 2);
            laguer(a, m, xr, xi, EPS, true);
            roots.host(2 * j + 1) = xr;
            roots.host(2 * j + 2) = xi;
        }
    }

    // Insertion sort by real part
    for (int j = 1; j < m; j++) {
        double kr = roots.host(2 * j + 1);
        double ki = roots.host(2 * j + 2);
        int i;
        for (i = j - 1; i >= 0; i--) {
            if (roots.host(2 * i + 1) <= kr) break;
            roots.host(2 * (i + 1) + 1) = roots.host(2 * i + 1);
            roots.host(2 * (i + 1) + 2) = roots.host(2 * i + 2);
        }
        roots.host(2 * (i + 1) + 1) = kr;
        roots.host(2 * (i + 1) + 2) = ki;
    }
}
