#include <cstdio>
#include <cmath>
#include <matar.h>
#include "ran3.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NSAMP = 10000;
        int idum = -1;
        Ran3State state;

        double sum = 0.0, sum2 = 0.0;
        double rmin = 1.0, rmax = 0.0;

        for (int i = 0; i < NSAMP; i++) {
            double x = ran3(idum, state);
            sum  += x;
            sum2 += x * x;
            if (x < rmin) rmin = x;
            if (x > rmax) rmax = x;
        }

        double mean = sum / NSAMP;
        double var  = sum2 / NSAMP - mean * mean;

        std::printf("ran3 — Knuth subtractive method\n");
        std::printf("  Samples : %d\n", NSAMP);
        std::printf("  Range   : [%.6f, %.6f]\n", rmin, rmax);
        std::printf("  Mean    : %.6f  (expected ~0.5)\n", mean);
        std::printf("  Variance: %.6f  (expected ~%.6f)\n", var, 1.0 / 12.0);

        bool pass = (rmin >= 0.0) && (rmax < 1.0) &&
                    (std::fabs(mean - 0.5) < 0.02) &&
                    (std::fabs(var - 1.0 / 12.0) < 0.005);
        std::printf("  Test    : %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
