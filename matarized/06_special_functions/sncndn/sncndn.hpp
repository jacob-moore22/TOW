#pragma once
#include <cmath>
#include <matar.h>

// Jacobian elliptic functions sn, cn, dn via arithmetic-geometric mean.
// emmc = 1 - m where m is the parameter (square of the modulus k).
KOKKOS_INLINE_FUNCTION
void sncndn(double uu, double emmc, double& sn, double& cn, double& dn)
{
    constexpr double CA = 0.0003;
    double em[13], en[13];
    double emc = emmc;
    double u = uu;

    if (emc != 0.0) {
        bool bo = (emc < 0.0);
        double d = 0.0;
        if (bo) {
            d = 1.0 - emc;
            emc = -emc / d;
            d = sqrt(d);
            u = d * u;
        }
        double a = 1.0;
        dn = 1.0;
        int l = 0;
        double c = 0.0;
        for (int i = 0; i < 13; i++) {
            l = i;
            em[i] = a;
            emc = sqrt(emc);
            en[i] = emc;
            c = 0.5 * (a + emc);
            if (fabs(a - emc) <= CA * a) break;
            emc = a * emc;
            a = c;
        }
        u = c * u;
        sn = sin(u);
        cn = cos(u);
        if (sn != 0.0) {
            a = cn / sn;
            c = a * c;
            for (int ii = l; ii >= 0; ii--) {
                double b = em[ii];
                a = c * a;
                c = dn * c;
                dn = (en[ii] + a) / (b + a);
                a = c / b;
            }
            a = 1.0 / sqrt(c * c + 1.0);
            if (sn < 0.0) {
                sn = -a;
            } else {
                sn = a;
            }
            cn = c * sn;
        }
        if (bo) {
            a = dn;
            dn = cn;
            cn = a;
            sn = sn / d;
        }
    } else {
        cn = 1.0 / cosh(u);
        dn = cn;
        sn = tanh(u);
    }
}
