#pragma once
#include <cstdio>
#include <cmath>
#include <matar.h>

using namespace mtr;

// Forward declarations of date dependencies
KOKKOS_INLINE_FUNCTION int  julday(int mm, int id, int iyyy);
KOKKOS_INLINE_FUNCTION void caldat(int julian, int& mm, int& id, int& iyyy);
KOKKOS_INLINE_FUNCTION void flmoon(int n, int nph, int& jd, double& frac);

// Friday the 13th full moon finder (Numerical Recipes BADLUK).
// Searches [iybeg, iyend] for Fridays the 13th that coincide with a full moon.
// timzon: timezone offset in fractional days (e.g. -5/24 for EST).
inline void badluk(int iybeg, int iyend, double timzon)
{
    std::printf("Full moons on Friday the 13th from %d to %d\n", iybeg, iyend);

    for (int iyyy = iybeg; iyyy <= iyend; iyyy++) {
        for (int im = 1; im <= 12; im++) {
            int jday = julday(im, 13, iyyy);
            int idwk = (jday + 1) % 7;

            if (idwk == 5) {
                int n = static_cast<int>(12.37 * (iyyy - 1900 + (im - 0.5) / 12.0));
                int icon = 0;

                for (;;) {
                    int jd;
                    double frac;
                    flmoon(n, 2, jd, frac);

                    int ifrac = static_cast<int>(24.0 * (frac + timzon) + 0.5);
                    if (ifrac < 0) {
                        jd = jd - 1;
                        ifrac = ifrac + 24;
                    }
                    if (ifrac > 12) {
                        jd = jd + 1;
                        ifrac = ifrac - 12;
                    } else {
                        ifrac = ifrac + 12;
                    }

                    if (jd == jday) {
                        std::printf("\n %d/%d/%d\n", im, 13, iyyy);
                        std::printf(" Full moon %d hrs after midnight (EST).\n", ifrac);
                        break;
                    } else {
                        int ic = (jday > jd) ? 1 : -1;
                        if (ic == -icon) break;
                        icon = ic;
                        n = n + ic;
                    }
                }
            }
        }
    }
}
