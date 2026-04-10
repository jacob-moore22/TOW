#include <cstdio>
#include <cmath>
#include <matar.h>
#include "flmoon.hpp"
#include "../caldat/caldat.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Full moon dates (first 12 full moons of 2000):\n");
        std::printf("%6s %12s %12s %10s\n", "N", "JD", "Date", "frac");

        // Approximate lunation number for Jan 2000
        int n_start = static_cast<int>(12.37 * (2000 - 1900));
        bool any_fail = false;

        for (int i = 0; i < 12; i++) {
            int jd;
            double frac;
            flmoon(n_start + i, 2, jd, frac);

            int mm, id, iyyy;
            caldat(jd, mm, id, iyyy);

            std::printf("%6d %12d %4d/%02d/%02d %10.4f\n",
                        n_start + i, jd, iyyy, mm, id, frac);

            if (iyyy < 1999 || iyyy > 2001) any_fail = true;
        }

        // Specific known full moon: Jan 21, 2000 (JD 2451564)
        int jd_check;
        double frac_check;
        flmoon(n_start, 2, jd_check, frac_check);
        int mm, id, iyyy;
        caldat(jd_check, mm, id, iyyy);
        std::printf("\nFirst full moon ~2000: %d/%d/%d (JD %d)\n",
                    mm, id, iyyy, jd_check);

        std::printf("Test %s\n", !any_fail ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
