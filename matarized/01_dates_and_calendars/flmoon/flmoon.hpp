#pragma once
#include <cmath>
#include <matar.h>

using namespace mtr;

// Full moon calculator (Numerical Recipes FLMOON).
// n: lunation number (0 = first new moon of 1900)
// nph: phase (0=new, 1=first quarter, 2=full, 3=last quarter)
// Outputs: jd (Julian Day of the phase), frac (fractional day past midnight).
KOKKOS_INLINE_FUNCTION
void flmoon(int n, int nph, int& jd, double& frac)
{
    constexpr double RAD = 0.017453293;

    double c  = n + nph / 4.0;
    double t  = c / 1236.85;
    double t2 = t * t;
    double as = 359.2242 + 29.105356 * c;
    double am = 306.0253 + 385.816918 * c + 0.010730 * t2;

    jd = 2415020 + 28 * n + 7 * nph;

    double xtra = 0.75933 + 1.53058868 * c + (1.178e-4 - 1.55e-7 * t) * t2;

    if (nph == 0 || nph == 2) {
        xtra += (0.1734 - 3.93e-4 * t) * sin(RAD * as)
              - 0.4068 * sin(RAD * am);
    } else if (nph == 1 || nph == 3) {
        xtra += (0.1721 - 4.0e-4 * t) * sin(RAD * as)
              - 0.6280 * sin(RAD * am);
    }

    int i;
    if (xtra >= 0.0)
        i = static_cast<int>(xtra);
    else
        i = static_cast<int>(xtra - 1.0);

    jd   = jd + i;
    frac = xtra - i;
}
