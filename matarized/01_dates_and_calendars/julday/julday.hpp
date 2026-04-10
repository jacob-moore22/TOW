#pragma once
#include <matar.h>

using namespace mtr;

// Julian Day Number from calendar date (Numerical Recipes JULDAY).
// mm: month (1-12), id: day, iyyy: year (negative for BC, no year 0).
// Returns -1 on error (year zero).
KOKKOS_INLINE_FUNCTION
int julday(int mm, int id, int iyyy)
{
    const int IGREG = 15 + 31 * (10 + 12 * 1582);

    if (iyyy == 0) return -1;
    if (iyyy < 0) iyyy = iyyy + 1;

    int jy, jm;
    if (mm > 2) {
        jy = iyyy;
        jm = mm + 1;
    } else {
        jy = iyyy - 1;
        jm = mm + 13;
    }

    int jday = static_cast<int>(365.25 * jy)
             + static_cast<int>(30.6001 * jm)
             + id + 1720995;

    if (id + 31 * (mm + 12 * iyyy) >= IGREG) {
        int ja = static_cast<int>(0.01 * jy);
        jday = jday + 2 - ja + static_cast<int>(0.25 * ja);
    }
    return jday;
}
