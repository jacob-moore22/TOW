#include <cstdio>
#include <cmath>
#include <matar.h>

using namespace mtr;

#include "../bksub/bksub.hpp"
#include "../pinvs/pinvs.hpp"
#include "../red/red.hpp"
#include "../solvde/solvde.hpp"
#include "../difeq/difeq.hpp"
#include "sfroid.hpp"

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        struct TestCase { int mm; int n; double c2; double expected; };
        TestCase cases[] = {
            {0, 2, 0.0, 6.0},
            {0, 4, 0.0, 20.0},
            {1, 2, 0.0, 6.0},
        };
        int ncases = 3;

        std::printf("Spheroidal wave eigenvalues (c2=0 should match l*(l+1)):\n");
        std::printf("%6s %6s %8s %16s %16s %12s\n",
                    "m", "n", "c^2", "lambda", "expected", "error");

        bool all_pass = true;
        for (int t = 0; t < ncases; t++) {
            double lambda = sfroid(cases[t].mm, cases[t].n, cases[t].c2);
            double error  = std::fabs(lambda - cases[t].expected);
            std::printf("%6d %6d %8.1f %16.10f %16.1f %12.2e\n",
                        cases[t].mm, cases[t].n, cases[t].c2,
                        lambda, cases[t].expected, error);
            if (error > 1e-3) all_pass = false;
        }

        std::printf("\nTest %s\n", all_pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
