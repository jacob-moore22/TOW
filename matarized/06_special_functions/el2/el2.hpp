#pragma once
#include <cmath>
#include <matar.h>

// General elliptic integral of the second kind.
// el2(x, kc, a, b) = integral from 0 to atan(x) of
//   (a + b*tan^2(phi)) / sqrt((1+tan^2(phi))*(1+kc^2*tan^2(phi))) dphi
KOKKOS_INLINE_FUNCTION
double el2(double x, double qqc, double aa, double bb)
{
    constexpr double PI = 3.14159265;
    constexpr double CA = 0.0003;
    constexpr double CB = 1.0e-9;

    if (x == 0.0) {
        return 0.0;
    } else if (qqc != 0.0) {
        double qc  = qqc;
        double a   = aa;
        double b   = bb;
        double c   = x * x;
        double d   = 1.0 + c;
        double p   = sqrt((1.0 + qc * qc * c) / d);
        d = x / d;
        c = d / (2.0 * p);
        double z   = a - b;
        double eye = a;
        a = 0.5 * (b + a);
        double y   = fabs(1.0 / x);
        double f   = 0.0;
        int    l   = 0;
        double em  = 1.0;
        qc = fabs(qc);

        for (;;) {
            b = eye * qc + b;
            double e = em * qc;
            double g = e / p;
            d = f * g + d;
            f = c;
            eye = a;
            p = g + p;
            c = 0.5 * (d / p + c);
            g = em;
            em = qc + em;
            a = 0.5 * (b / em + a);
            y = -e / y + y;
            if (y == 0.0) y = sqrt(e) * CB;
            if (fabs(g - qc) > CA * g) {
                qc = sqrt(e) * 2.0;
                l = l + l;
                if (y < 0.0) l = l + 1;
            } else {
                break;
            }
        }

        if (y < 0.0) l = l + 1;
        double result = (atan(em / y) + PI * l) * a / em;
        if (x < 0.0) result = -result;
        return result + c * z;
    } else {
        return 0.0;
    }
}
