#include <cstdio>
#include <cmath>
#include <matar.h>
#include "gamdev.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NSAMP = 10000;
        int idum = -1;
        Ran1State state;

        int test_orders[] = {1, 3, 5, 10};
        int norders = 4;

        std::printf("gamdev — Gamma deviates of integer order\n");
        std::printf("%8s %12s %12s %12s %12s\n",
                    "Order", "Mean", "Exp.Mean", "Var", "Exp.Var");

        bool all_pass = true;
        for (int t = 0; t < norders; t++) {
            int ia = test_orders[t];
            double sum = 0.0, sum2 = 0.0;
            for (int i = 0; i < NSAMP; i++) {
                double x = gamdev(ia, idum, state);
                sum  += x;
                sum2 += x * x;
            }
            double mean = sum / NSAMP;
            double var  = sum2 / NSAMP - mean * mean;
            double exp_mean = static_cast<double>(ia);
            double exp_var  = static_cast<double>(ia);

            std::printf("%8d %12.4f %12.4f %12.4f %12.4f\n",
                        ia, mean, exp_mean, var, exp_var);

            if (std::fabs(mean - exp_mean) > 0.15 * exp_mean) all_pass = false;
            if (std::fabs(var  - exp_var)  > 0.25 * exp_var)  all_pass = false;
        }
        std::printf("  Test : %s\n", all_pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
