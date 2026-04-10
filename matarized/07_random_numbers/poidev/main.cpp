#include <cstdio>
#include <cmath>
#include <matar.h>
#include "poidev.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NSAMP = 10000;
        int idum = -1;
        Ran1State rstate;
        PoidevState pstate;

        double test_xm[] = {1.0, 5.0, 20.0, 100.0};
        int ntest = 4;

        std::printf("poidev — Poisson deviates\n");
        std::printf("%10s %12s %12s %12s %12s\n",
                    "xm", "Mean", "Exp.Mean", "Var", "Exp.Var");

        bool all_pass = true;
        for (int t = 0; t < ntest; t++) {
            double xm = test_xm[t];
            pstate = PoidevState();

            double sum = 0.0, sum2 = 0.0;
            for (int i = 0; i < NSAMP; i++) {
                double x = poidev(xm, idum, rstate, pstate);
                sum  += x;
                sum2 += x * x;
            }
            double mean = sum / NSAMP;
            double var  = sum2 / NSAMP - mean * mean;

            std::printf("%10.1f %12.4f %12.4f %12.4f %12.4f\n",
                        xm, mean, xm, var, xm);

            if (std::fabs(mean - xm) > 0.15 * xm + 0.5) all_pass = false;
            if (std::fabs(var  - xm) > 0.25 * xm + 1.0) all_pass = false;
        }
        std::printf("  Test : %s\n", all_pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
