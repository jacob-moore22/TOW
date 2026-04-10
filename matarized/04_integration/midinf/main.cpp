#include <cstdio>
#include <cmath>
#include <matar.h>
#include "midinf.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Midpoint rule for semi-infinite intervals (midinf)\n\n");

        // Test: integral of exp(-x) from 1 to infinity = 1/e
        auto f_exp_neg = [](double x) { return std::exp(-x); };
        double expected = std::exp(-1.0);

        std::printf("exp(-x) on [1, inf), exact = 1/e = %.12f\n\n", expected);
        std::printf("%4s %18s %14s\n", "n", "Estimate", "Error");

        double s;
        int it = 0;
        for (int n = 1; n <= 15; n++) {
            midinf(f_exp_neg, 1.0, 1.0e30, s, n, it);
            std::printf("%4d %18.12f %14.2e\n", n, s, std::fabs(s - expected));
        }

        double err1 = std::fabs(s - expected);
        std::printf("\nFinal: result=%.12f  exact=%.12f  err=%.2e\n",
                    s, expected, err1);

        // Test: integral of 1/x^2 from 1 to infinity = 1
        auto f_invx2 = [](double x) { return 1.0 / (x * x); };
        double expected2 = 1.0;

        double s2;
        int it2 = 0;
        for (int n = 1; n <= 15; n++)
            midinf(f_invx2, 1.0, 1.0e30, s2, n, it2);

        double err2 = std::fabs(s2 - expected2);
        std::printf("\n1/x^2 on [1, inf), exact = 1:\n");
        std::printf("  Result   = %.12f\n", s2);
        std::printf("  Expected = %.12f\n", expected2);
        std::printf("  Error    = %.2e\n", err2);

        bool pass = (err1 < 1e-4) && (err2 < 1e-4);
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
