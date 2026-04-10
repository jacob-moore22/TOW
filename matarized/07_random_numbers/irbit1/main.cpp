#include <cstdio>
#include <cmath>
#include <matar.h>
#include "irbit1.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NSAMP = 10000;
        int iseed = 12345;

        int count1 = 0;
        for (int i = 0; i < NSAMP; i++) {
            count1 += irbit1(iseed);
        }

        double frac = static_cast<double>(count1) / NSAMP;

        std::printf("irbit1 — LFSR random bit generator (taps 1,2,5,18)\n");
        std::printf("  Samples : %d\n", NSAMP);
        std::printf("  Ones    : %d  (%.4f)\n", count1, frac);
        std::printf("  Zeros   : %d  (%.4f)\n", NSAMP - count1, 1.0 - frac);

        bool pass = std::fabs(frac - 0.5) < 0.02;
        std::printf("  Test    : %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
