#pragma once
#include <matar.h>

using namespace mtr;

// Calendar date from Julian Day Number (Numerical Recipes CALDAT).
// julian: Julian Day Number
// Outputs: mm (month), id (day), iyyy (year).
KOKKOS_INLINE_FUNCTION
void caldat(int julian, int& mm, int& id, int& iyyy)
{
    const int IGREG = 2299161;

    int ja;
    if (julian >= IGREG) {
        int jalpha = static_cast<int>(((julian - 1867216) - 0.25) / 36524.25);
        ja = julian + 1 + jalpha - static_cast<int>(0.25 * jalpha);
    } else {
        ja = julian;
    }

    int jb = ja + 1524;
    int jc = static_cast<int>(6680.0 + ((jb - 2439870) - 122.1) / 365.25);
    int jd = 365 * jc + static_cast<int>(0.25 * jc);
    int je = static_cast<int>((jb - jd) / 30.6001);

    id = jb - jd - static_cast<int>(30.6001 * je);
    mm = je - 1;
    if (mm > 12) mm = mm - 12;
    iyyy = jc - 4715;
    if (mm > 2) iyyy = iyyy - 1;
    if (iyyy <= 0) iyyy = iyyy - 1;
}
