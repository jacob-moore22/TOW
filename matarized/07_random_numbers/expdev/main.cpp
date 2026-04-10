#include <cstdio>
#include <cmath>
#include <matar.h>
#include "expdev.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NSAMP = 10000;
        int idum = -1;
        Ran1State state;

        double sum = 0.0, sum2 = 0.0;
        double rmin = 1e30, rmax = -1e30;

        for (int i = 0; i < NSAMP; i++) {
            double x = expdev(idum, state);
            sum  += x;
            sum2 += x * x;
            if (x < rmin) rmin = x;
            if (x > rmax) rmax = x;
        }

        double mean = sum / NSAMP;
        double var  = sum2 / NSAMP - mean * mean;

        std::printf("expdev — exponential deviates (unit mean)\n");
        std::printf("  Samples : %d\n", NSAMP);
        std::printf("  Range   : [%.6f, %.6f]\n", rmin, rmax);
        std::printf("  Mean    : %.6f  (expected ~1.0)\n", mean);
        std::printf("  Variance: %.6f  (expected ~1.0)\n", var);

        bool pass = (rmin >= 0.0) &&
                    (std::fabs(mean - 1.0) < 0.05) &&
                    (std::fabs(var - 1.0) < 0.1);
        std::printf("  Test    : %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
