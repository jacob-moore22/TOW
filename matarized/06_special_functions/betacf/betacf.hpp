#pragma once
#include <cmath>
#include <matar.h>

// Continued-fraction evaluation for the incomplete beta function I_x(a,b).
// Used by betai; not called directly.
KOKKOS_INLINE_FUNCTION
double betacf(double a, double b, double x)
{
    constexpr int ITMAX = 100;
    constexpr double EPS = 3.0e-7;

    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double bz  = 1.0 - qab * x / qap;

    double am = 1.0;
    double bm = 1.0;
    double az = 1.0;

    for (int m = 1; m <= ITMAX; m++) {
        double em  = double(m);
        double tem = em + em;

        double d  = em * (b - em) * x / ((qam + tem) * (a + tem));
        double ap = az + d * am;
        double bp = bz + d * bm;

        d = -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem));
        double app = ap + d * az;
        double bpp = bp + d * bz;

        double aold = az;
        am = ap / bpp;
        bm = bp / bpp;
        az = app / bpp;
        bz = 1.0;

        if (fabs(az - aold) < EPS * fabs(az)) {
            return az;
        }
    }
    return az;
}
