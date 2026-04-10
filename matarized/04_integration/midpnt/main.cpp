#include <cstdio>
#include <cmath>
#include <matar.h>
#include "midpnt.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const double pi = 3.14159265358979323846;

        std::printf("Midpoint rule refinement tests\n\n");

        // Test 1: integral of sin(x) from 0 to pi = 2
        auto f_sin = [](double x) { return std::sin(x); };
        double s1;
        int it1 = 0;
        std::printf("sin(x) on [0, pi], exact = 2.0:\n");
        std::printf("%4s %18s %14s\n", "n", "Estimate", "Error");
        for (int n = 1; n <= 12; n++) {
            midpnt(f_sin, 0.0, pi, s1, n, it1);
            std::printf("%4d %18.12f %14.2e\n", n, s1, std::fabs(s1 - 2.0));
        }

        // Test 2: integral of x^2 from 0 to 1 = 1/3
        auto f_x2 = [](double x) { return x * x; };
        double s2;
        int it2 = 0;
        for (int n = 1; n <= 12; n++)
            midpnt(f_x2, 0.0, 1.0, s2, n, it2);

        std::printf("\nx^2 on [0,1], exact = 1/3:\n");
        std::printf("  Result   = %.12f\n", s2);
        std::printf("  Expected = %.12f\n", 1.0 / 3.0);
        std::printf("  Error    = %.2e\n", std::fabs(s2 - 1.0 / 3.0));

        // Test 3: integral of exp(x) from 0 to 1 = e-1
        auto f_exp = [](double x) { return std::exp(x); };
        double s3;
        int it3 = 0;
        for (int n = 1; n <= 12; n++)
            midpnt(f_exp, 0.0, 1.0, s3, n, it3);

        double e_val = std::exp(1.0) - 1.0;
        std::printf("\nexp(x) on [0,1], exact = e-1:\n");
        std::printf("  Result   = %.12f\n", s3);
        std::printf("  Expected = %.12f\n", e_val);
        std::printf("  Error    = %.2e\n", std::fabs(s3 - e_val));

        bool pass = (std::fabs(s1 - 2.0) < 1e-4) &&
                    (std::fabs(s2 - 1.0 / 3.0) < 1e-4) &&
                    (std::fabs(s3 - e_val) < 1e-4);
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
