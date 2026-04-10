#include <cstdio>
#include <cmath>
#include <matar.h>
#include "gasdev.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NSAMP = 10000;
        int idum = -1;
        Ran1State rstate;
        GasdevState gstate;

        double sum = 0.0, sum2 = 0.0;

        for (int i = 0; i < NSAMP; i++) {
            double x = gasdev(idum, rstate, gstate);
            sum  += x;
            sum2 += x * x;
        }

        double mean = sum / NSAMP;
        double var  = sum2 / NSAMP - mean * mean;

        std::printf("gasdev — Gaussian deviates (Box-Muller)\n");
        std::printf("  Samples : %d\n", NSAMP);
        std::printf("  Mean    : %.6f  (expected ~0.0)\n", mean);
        std::printf("  Variance: %.6f  (expected ~1.0)\n", var);

        bool pass = (std::fabs(mean) < 0.05) &&
                    (std::fabs(var - 1.0) < 0.1);
        std::printf("  Test    : %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
