#pragma once
#include <cmath>
#include <complex>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// Laguerre's method for finding a root of a polynomial with complex
// coefficients. a[0..m] are the coefficients (a[0] = constant term,
// a[m] = leading coefficient). x is the initial guess, refined in-place.
// If polish is true, the convergence criterion uses EPS rounding only.
inline void laguer(const std::complex<double>* a, int m,
                   std::complex<double>& x, double eps, bool polish)
{
    constexpr double EPSS  = 6.0e-8;
    constexpr int    MAXIT = 100;

    double dxold = std::abs(x);

    for (int iter = 0; iter < MAXIT; iter++) {
        std::complex<double> b = a[m];
        double err = std::abs(b);
        std::complex<double> d(0.0, 0.0);
        std::complex<double> f(0.0, 0.0);
        double abx = std::abs(x);

        for (int j = m - 1; j >= 0; j--) {
            f = x * f + d;
            d = x * d + b;
            b = x * b + a[j];
            err = std::abs(b) + abx * err;
        }
        err *= EPSS;

        if (std::abs(b) <= err) return;

        std::complex<double> g = d / b;
        std::complex<double> g2 = g * g;
        std::complex<double> h = g2 - 2.0 * f / b;
        std::complex<double> sq = std::sqrt(
            static_cast<double>(m - 1) * (static_cast<double>(m) * h - g2));
        std::complex<double> gp = g + sq;
        std::complex<double> gm = g - sq;
        if (std::abs(gp) < std::abs(gm)) gp = gm;
        std::complex<double> dx = static_cast<double>(m) / gp;

        std::complex<double> x1 = x - dx;
        if (x.real() == x1.real() && x.imag() == x1.imag()) return;
        x = x1;

        double cdx = std::abs(dx);
        if (iter > 6 && cdx >= dxold) return;
        dxold = cdx;

        if (!polish) {
            if (std::abs(dx) <= eps * std::abs(x)) return;
        }
    }

    std::fprintf(stderr, "laguer: too many iterations\n");
}
