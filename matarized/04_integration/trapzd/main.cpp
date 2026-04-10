#include <cstdio>
#include <cmath>
#include <matar.h>
#include "trapzd.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const double pi       = 3.14159265358979323846;
        const double expected = 2.0; // integral of sin(x) from 0 to pi

        auto f_sin = [](double x) { return std::sin(x); };

        std::printf("Trapezoidal rule refinement for integral of sin(x) from 0 to pi\n");
        std::printf("Exact value = 2.0\n\n");
        std::printf("%4s %18s %14s\n", "n", "Estimate", "Error");

        double s;
        int it = 0;
        double max_err = 0.0;

        for (int n = 1; n <= 15; n++) {
            trapzd(f_sin, 0.0, pi, s, n, it);
            double err = std::fabs(s - expected);
            if (err > max_err) max_err = err;
            std::printf("%4d %18.12f %14.2e\n", n, s, err);
        }

        // Also test x^2 on [0,1], exact = 1/3
        double s2;
        int it2 = 0;
        auto f_x2 = [](double x) { return x * x; };
        for (int n = 1; n <= 15; n++)
            trapzd(f_x2, 0.0, 1.0, s2, n, it2);

        double err2 = std::fabs(s2 - 1.0 / 3.0);
        std::printf("\nIntegral of x^2 from 0 to 1 (after 15 refinements):\n");
        std::printf("  Result   = %.12f\n", s2);
        std::printf("  Expected = %.12f\n", 1.0 / 3.0);
        std::printf("  Error    = %.2e\n", err2);

        bool pass = (std::fabs(s - expected) < 1e-4) && (err2 < 1e-4);
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
