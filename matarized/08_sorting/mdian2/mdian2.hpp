#pragma once
#include <cmath>
#include <algorithm>
#include <matar.h>

using namespace mtr;

// Median via selection. Returns the median of x(1..n) without fully sorting
// the array, using an iterative estimation scheme.
inline void mdian2(DFMatrixKokkos<double>& x, int n, double& xmed)
{
    constexpr double BIG  = 1.0e30;
    constexpr double AFAC = 1.5;
    constexpr double AMP  = 1.5;

    double a   = 0.5 * (x.host(1) + x.host(n));
    double eps = std::fabs(x.host(n) - x.host(1));
    double ap  =  BIG;
    double am  = -BIG;

    while (true) {
        double sum  = 0.0;
        double sumx = 0.0;
        int    np   = 0;
        int    nm   = 0;
        double xp   =  BIG;
        double xm   = -BIG;

        for (int j = 1; j <= n; j++) {
            double xx = x.host(j);
            if (xx != a) {
                if (xx > a) {
                    np++;
                    if (xx < xp) xp = xx;
                } else if (xx < a) {
                    nm++;
                    if (xx > xm) xm = xx;
                }
                double dum = 1.0 / (eps + std::fabs(xx - a));
                sum  += dum;
                sumx += xx * dum;
            }
        }

        if (np - nm >= 2) {
            am = a;
            double aa = xp + std::max(0.0, sumx / sum - a) * AMP;
            if (aa > ap) aa = 0.5 * (a + ap);
            eps = AFAC * std::fabs(aa - a);
            a   = aa;
        } else if (nm - np >= 2) {
            ap = a;
            double aa = xm + std::min(0.0, sumx / sum - a) * AMP;
            if (aa < am) aa = 0.5 * (a + am);
            eps = AFAC * std::fabs(aa - a);
            a   = aa;
        } else {
            if (n % 2 == 0) {
                if (np == nm) {
                    xmed = 0.5 * (xp + xm);
                } else if (np > nm) {
                    xmed = 0.5 * (a + xp);
                } else {
                    xmed = 0.5 * (xm + a);
                }
            } else {
                if (np == nm) {
                    xmed = a;
                } else if (np > nm) {
                    xmed = xp;
                } else {
                    xmed = xm;
                }
            }
            return;
        }
    }
}
