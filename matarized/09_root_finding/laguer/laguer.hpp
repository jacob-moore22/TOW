// Laguerre's method for finding a root of a polynomial with complex
// coefficients (Numerical Recipes LAGUER).
//
// a:      interleaved real/imag coefficient array, DFMatrixKokkos<double>(2*(m+1)).
//         Coefficient j (0-based) is stored as:
//           real part at a(2*j + 1), imag part at a(2*j + 2)   [1-based indexing]
// m:      polynomial degree
// xr,xi:  real and imaginary parts of root -- initial guess in, refined result out
// eps:    convergence tolerance on |dx|/|x|
// polish: if true, iterate to full machine precision (ignore eps)
//
// Both loops are inherently serial:
//   - Outer iteration: each Laguerre step depends on the previous x
//   - Inner Horner evaluation: b(j) depends on b(j+1) (loop-carried dependency)
// Parallelism for polynomial root-finding lives at the caller level
// (e.g., zroots finding all M roots independently).

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

// GPU-portable complex number (plain struct, no std:: dependency)
struct PortableComplex {
    double re, im;
    KOKKOS_INLINE_FUNCTION PortableComplex() : re(0.0), im(0.0) {}
    KOKKOS_INLINE_FUNCTION PortableComplex(double r, double i = 0.0) : re(r), im(i) {}
    KOKKOS_INLINE_FUNCTION PortableComplex operator+(PortableComplex o) const {
        return {re + o.re, im + o.im};
    }
    KOKKOS_INLINE_FUNCTION PortableComplex operator-(PortableComplex o) const {
        return {re - o.re, im - o.im};
    }
    KOKKOS_INLINE_FUNCTION PortableComplex operator*(PortableComplex o) const {
        return {re * o.re - im * o.im, re * o.im + im * o.re};
    }
    KOKKOS_INLINE_FUNCTION PortableComplex operator/(PortableComplex o) const {
        double d = o.re * o.re + o.im * o.im;
        return {(re * o.re + im * o.im) / d, (im * o.re - re * o.im) / d};
    }
    KOKKOS_INLINE_FUNCTION PortableComplex& operator-=(PortableComplex o) {
        re -= o.re; im -= o.im; return *this;
    }
};

KOKKOS_INLINE_FUNCTION
double pc_abs(PortableComplex z) { return sqrt(z.re * z.re + z.im * z.im); }

KOKKOS_INLINE_FUNCTION
PortableComplex pc_sqrt(PortableComplex z) {
    double r = pc_abs(z);
    if (r == 0.0) return {0.0, 0.0};
    if (z.re >= 0.0) {
        double s = sqrt(0.5 * (r + z.re));
        return {s, z.im / (2.0 * s)};
    }
    double s = sqrt(0.5 * (r - z.re));
    if (z.im < 0.0) s = -s;
    return {z.im / (2.0 * s), s};
}

KOKKOS_INLINE_FUNCTION
PortableComplex pc_conj(PortableComplex z) { return {z.re, -z.im}; }

inline void laguer(DFMatrixKokkos<double>& a, int m,
                   double& xr, double& xi,
                   double eps, bool polish)
{
    constexpr double EPSS  = 6.0e-8;
    constexpr int    MAXIT = 100;
    constexpr int    MAX_DEGREE = 128;

    PortableComplex coef[MAX_DEGREE + 1];
    for (int j = 0; j <= m; j++)
        coef[j] = {a.host(2 * j + 1), a.host(2 * j + 2)};

    PortableComplex x(xr, xi);
    double dxold = pc_abs(x);

    for (int iter = 0; iter < MAXIT; iter++) {
        PortableComplex b = coef[m];
        double err = pc_abs(b);
        PortableComplex d(0.0, 0.0), f(0.0, 0.0);
        double abx = pc_abs(x);

        for (int j = m - 1; j >= 0; j--) {
            f = x * f + d;
            d = x * d + b;
            b = x * b + coef[j];
            err = pc_abs(b) + abx * err;
        }
        err *= EPSS;

        if (pc_abs(b) <= err) { xr = x.re; xi = x.im; return; }

        PortableComplex g  = d / b;
        PortableComplex g2 = g * g;
        PortableComplex h  = g2 - PortableComplex(2.0) * f / b;
        PortableComplex sq = pc_sqrt(
            PortableComplex(static_cast<double>(m - 1)) *
            (PortableComplex(static_cast<double>(m)) * h - g2));
        PortableComplex gp = g + sq;
        PortableComplex gm = g - sq;
        if (pc_abs(gp) < pc_abs(gm)) gp = gm;
        PortableComplex dx = PortableComplex(static_cast<double>(m)) / gp;

        PortableComplex x1 = x - dx;
        if (x.re == x1.re && x.im == x1.im) { xr = x.re; xi = x.im; return; }
        x = x1;

        double cdx = pc_abs(dx);
        if (iter > 6 && cdx >= dxold) { xr = x.re; xi = x.im; return; }
        dxold = cdx;

        if (!polish) {
            if (pc_abs(dx) <= eps * pc_abs(x)) { xr = x.re; xi = x.im; return; }
        }
    }

    xr = x.re; xi = x.im;
    std::fprintf(stderr, "laguer: too many iterations\n");
}
