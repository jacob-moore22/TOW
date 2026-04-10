#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bnldev.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NSAMP = 10000;
        int idum = -1;
        Ran1State rstate;
        BnldevState bstate;

        struct TestCase { double pp; int n; };
        TestCase cases[] = {{0.3, 10}, {0.5, 20}, {0.7, 50}, {0.1, 100}};
        int ncases = 4;

        std::printf("bnldev — Binomial deviates\n");
        std::printf("%6s %6s %12s %12s %12s %12s\n",
                    "pp", "n", "Mean", "Exp.Mean", "Var", "Exp.Var");

        bool all_pass = true;
        for (int t = 0; t < ncases; t++) {
            double pp = cases[t].pp;
            int    n  = cases[t].n;
            bstate = BnldevState();

            double sum = 0.0, sum2 = 0.0;
            for (int i = 0; i < NSAMP; i++) {
                double x = bnldev(pp, n, idum, rstate, bstate);
                sum  += x;
                sum2 += x * x;
            }
            double mean     = sum / NSAMP;
            double var      = sum2 / NSAMP - mean * mean;
            double exp_mean = n * pp;
            double exp_var  = n * pp * (1.0 - pp);

            std::printf("%6.2f %6d %12.4f %12.4f %12.4f %12.4f\n",
                        pp, n, mean, exp_mean, var, exp_var);

            if (std::fabs(mean - exp_mean) > 0.15 * exp_mean + 0.5) all_pass = false;
            if (std::fabs(var  - exp_var)  > 0.30 * exp_var  + 0.5) all_pass = false;
        }
        std::printf("  Test : %s\n", all_pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
