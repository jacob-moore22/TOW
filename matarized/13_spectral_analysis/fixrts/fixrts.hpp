// Fix unstable linear-prediction roots by reflecting them inside
// the unit circle.  Self-contained: includes Laguerre root-finder
// and polynomial deflation (host-side, since npoles is small).
//
// d:      LP coefficients, npoles elements (1-based DFMatrixKokkos).
//         Modified in-place so all roots of the characteristic
//         polynomial lie inside the unit circle.
// npoles: number of LP poles.
//
// Host-side computation: npoles is typically O(10-100), so
// Laguerre iteration + deflation is run on the host.

#pragma once
#include <cmath>
#include <cstdio>
#include <matar.h>

using namespace mtr;

namespace fixrts_detail {

// GPU-portable complex number (same as laguer.hpp PortableComplex)
struct FC {
    double re, im;
    FC() : re(0.0), im(0.0) {}
    FC(double r, double i = 0.0) : re(r), im(i) {}
    FC operator+(FC o) const { return {re + o.re, im + o.im}; }
    FC operator-(FC o) const { return {re - o.re, im - o.im}; }
    FC operator*(FC o) const { return {re*o.re - im*o.im, re*o.im + im*o.re}; }
    FC operator/(FC o) const {
        double d = o.re*o.re + o.im*o.im;
        return {(re*o.re + im*o.im)/d, (im*o.re - re*o.im)/d};
    }
    FC& operator-=(FC o) { re -= o.re; im -= o.im; return *this; }
};

inline double fc_abs(FC z) { return sqrt(z.re*z.re + z.im*z.im); }

inline FC fc_sqrt(FC z) {
    double r = fc_abs(z);
    if (r == 0.0) return {0.0, 0.0};
    if (z.re >= 0.0) {
        double s = sqrt(0.5 * (r + z.re));
        return {s, z.im / (2.0 * s)};
    }
    double s = sqrt(0.5 * (r - z.re));
    if (z.im < 0.0) s = -s;
    return {z.im / (2.0 * s), s};
}

inline FC fc_conj(FC z) { return {z.re, -z.im}; }

inline void laguer(FC* a, int m, FC& x, double eps, bool /*polish*/)
{
    const int MAXIT = 80;
    for (int iter = 1; iter <= MAXIT; iter++) {
        FC b = a[m];
        double err = fc_abs(b);
        FC d(0.0, 0.0), f(0.0, 0.0);
        double abx = fc_abs(x);
        for (int j = m - 1; j >= 0; j--) {
            f = x * f + d;
            d = x * d + b;
            b = x * b + a[j];
            err = fc_abs(b) + abx * err;
        }
        err *= eps;
        if (fc_abs(b) <= err) return;
        FC g  = d / b;
        FC g2 = g * g;
        FC h  = g2 - FC(2.0) * f / b;
        FC sq = fc_sqrt(FC(static_cast<double>(m - 1)) *
                        (FC(static_cast<double>(m)) * h - g2));
        FC gp = g + sq;
        FC gm = g - sq;
        if (fc_abs(gp) < fc_abs(gm)) gp = gm;
        FC dx;
        if (fc_abs(gp) > 0.0)
            dx = FC(static_cast<double>(m)) / gp;
        else
            dx = FC((1.0 + abx) * cos(static_cast<double>(iter)),
                    (1.0 + abx) * sin(static_cast<double>(iter)));
        x -= dx;
        if (fc_abs(dx) <= eps * fc_abs(x)) return;
    }
}

inline void zroots(FC* a, int m, FC* roots, bool polish)
{
    const double eps = 1.0e-6;
    const int MAXM = 102;
    FC ad[MAXM];

    for (int j = 0; j <= m; j++) ad[j] = a[j];

    for (int j = m; j >= 1; j--) {
        FC x(0.0, 0.0);
        laguer(ad, j, x, eps, false);
        if (fabs(x.im) <= 2.0 * eps * eps * fabs(x.re))
            x = FC(x.re, 0.0);
        roots[j - 1] = x;
        FC b = ad[j];
        for (int jj = j - 1; jj >= 0; jj--) {
            FC c = ad[jj];
            ad[jj] = b;
            b = x * b + c;
        }
    }

    if (polish) {
        for (int j = 0; j < m; j++)
            laguer(a, m, roots[j], eps, true);
    }

    for (int j = 1; j < m; j++) {
        FC key = roots[j];
        int i = j - 1;
        while (i >= 0 && roots[i].re > key.re) {
            roots[i + 1] = roots[i];
            i--;
        }
        roots[i + 1] = key;
    }
}

} // namespace fixrts_detail

inline void fixrts(DFMatrixKokkos<double>& d, int npoles)
{
    using fixrts_detail::FC;
    using fixrts_detail::fc_abs;
    using fixrts_detail::fc_conj;

    d.update_host();

    const int NPMAX = 101;
    FC a[NPMAX], roots[NPMAX];

    a[npoles] = FC(1.0, 0.0);
    for (int j = npoles - 1; j >= 0; j--)
        a[j] = FC(-d.host(npoles - j), 0.0);

    fixrts_detail::zroots(a, npoles, roots, true);

    for (int j = 0; j < npoles; j++) {
        if (fc_abs(roots[j]) > 1.0)
            roots[j] = FC(1.0) / fc_conj(roots[j]);
    }

    a[0] = FC(0.0) - roots[0];
    a[1] = FC(1.0, 0.0);
    for (int j = 1; j < npoles; j++) {
        a[j + 1] = FC(1.0, 0.0);
        for (int i = j; i >= 1; i--)
            a[i] = a[i - 1] - roots[j] * a[i];
        a[0] = FC(0.0) - roots[j] * a[0];
    }

    for (int j = 0; j < npoles; j++)
        d.host(npoles - j) = -a[j].re;

    d.update_device();
}
