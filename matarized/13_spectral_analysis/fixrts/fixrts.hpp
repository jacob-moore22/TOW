// Fix unstable linear-prediction roots by reflecting them inside
// the unit circle.  Self-contained: includes Laguerre root-finder
// and polynomial deflation (host-side, since npoles is small).
//
// d:      LP coefficients, npoles elements (1-based DFMatrixKokkos).
//         Modified in-place so all roots of the characteristic
//         polynomial lie inside the unit circle.
// npoles: number of LP poles.
//
// Host-side computation: npoles is typically O(10–100), so
// Laguerre iteration + deflation is run on the host.

#pragma once
#include <cmath>
#include <complex>
#include <algorithm>
#include <matar.h>

using namespace mtr;

namespace fixrts_detail {

inline void laguer(std::complex<double>* a, int m,
                   std::complex<double>& x, double eps, bool /*polish*/)
{
    const int MAXIT = 80;
    for (int iter = 1; iter <= MAXIT; iter++) {
        std::complex<double> b = a[m];
        double err = std::abs(b);
        std::complex<double> d(0.0, 0.0), f(0.0, 0.0);
        double abx = std::abs(x);
        for (int j = m - 1; j >= 0; j--) {
            f = x * f + d;
            d = x * d + b;
            b = x * b + a[j];
            err = std::abs(b) + abx * err;
        }
        err *= eps;
        if (std::abs(b) <= err) return;
        std::complex<double> g  = d / b;
        std::complex<double> g2 = g * g;
        std::complex<double> h  = g2 - 2.0 * f / b;
        std::complex<double> sq = std::sqrt(
            static_cast<double>(m - 1) * (static_cast<double>(m) * h - g2));
        std::complex<double> gp = g + sq;
        std::complex<double> gm = g - sq;
        if (std::abs(gp) < std::abs(gm)) gp = gm;
        std::complex<double> dx;
        if (std::abs(gp) > 0.0)
            dx = static_cast<double>(m) / gp;
        else
            dx = std::polar(1.0 + abx, static_cast<double>(iter));
        x -= dx;
        if (std::abs(dx) <= eps * std::abs(x)) return;
    }
}

inline void zroots(std::complex<double>* a, int m,
                   std::complex<double>* roots, bool polish)
{
    const double eps = 1.0e-6;
    const int MAXM = 102;
    std::complex<double> ad[MAXM];

    for (int j = 0; j <= m; j++) ad[j] = a[j];

    for (int j = m; j >= 1; j--) {
        std::complex<double> x(0.0, 0.0);
        laguer(ad, j, x, eps, false);
        if (std::abs(x.imag()) <= 2.0 * eps * eps * std::abs(x.real()))
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
            laguer(a, m, roots[j], eps, true);
    }

    // Insertion sort by real part
    for (int j = 1; j < m; j++) {
        std::complex<double> key = roots[j];
        int i = j - 1;
        while (i >= 0 && roots[i].real() > key.real()) {
            roots[i + 1] = roots[i];
            i--;
        }
        roots[i + 1] = key;
    }
}

} // namespace fixrts_detail

inline void fixrts(DFMatrixKokkos<double>& d, int npoles)
{
    d.update_host();

    const int NPMAX = 101;
    std::complex<double> a[NPMAX], roots[NPMAX];

    // Build characteristic polynomial from LP coefficients
    a[npoles] = std::complex<double>(1.0, 0.0);
    for (int j = npoles - 1; j >= 0; j--)
        a[j] = std::complex<double>(-d.host(npoles - j), 0.0);

    fixrts_detail::zroots(a, npoles, roots, true);

    // Reflect any root outside the unit circle
    for (int j = 0; j < npoles; j++) {
        if (std::abs(roots[j]) > 1.0)
            roots[j] = 1.0 / std::conj(roots[j]);
    }

    // Reconstruct polynomial from modified roots
    a[0] = -roots[0];
    a[1] = std::complex<double>(1.0, 0.0);
    for (int j = 1; j < npoles; j++) {
        a[j + 1] = std::complex<double>(1.0, 0.0);
        for (int i = j; i >= 1; i--)
            a[i] = a[i - 1] - roots[j] * a[i];
        a[0] = -roots[j] * a[0];
    }

    // Write stabilised LP coefficients back
    for (int j = 0; j < npoles; j++)
        d.host(npoles - j) = -a[j].real();

    d.update_device();
}
