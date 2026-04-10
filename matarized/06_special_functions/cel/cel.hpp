#pragma once
#include <cmath>
#include <matar.h>

// General complete elliptic integral via arithmetic-geometric mean.
// cel(kc, p, a, b) = integral from 0 to pi/2 of
//   (a*cos^2 + b*sin^2) / ((cos^2 + p*sin^2) * sqrt(cos^2 + kc^2*sin^2)) dtheta
KOKKOS_INLINE_FUNCTION
double cel(double qqc, double pp, double aa, double bb)
{
    constexpr double CA  = 0.0003;
    constexpr double PIO2 = 1.5707963268;

    if (qqc == 0.0) return 0.0;

    double qc = fabs(qqc);
    double a  = aa;
    double b  = bb;
    double p  = pp;
    double e  = qc;
    double em = 1.0;

    if (p > 0.0) {
        p = sqrt(p);
        b = b / p;
    } else {
        double f = qc * qc;
        double q = 1.0 - f;
        double g = 1.0 - p;
        f = f - p;
        q = q * (b - a * p);
        p = sqrt(f / g);
        a = (a - b) / g;
        b = -q / (g * g * p) + a * p;
    }

    for (;;) {
        double f = a;
        a = a + b / p;
        double g = e / p;
        b = b + f * g;
        b = b + b;
        p = g + p;
        g = em;
        em = qc + em;
        if (fabs(g - qc) > g * CA) {
            qc = sqrt(e);
            qc = qc + qc;
            e = qc * em;
        } else {
            break;
        }
    }

    return PIO2 * (b + a * em) / (em * (em + p));
}
