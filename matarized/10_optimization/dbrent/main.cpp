#include <cstdio>
#include <cmath>
#include <matar.h>
#include "dbrent.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Test 1: f(x) = (x-3)^2, f'(x) = 2(x-3), minimum at x = 3
        auto f  = [](double x) { return (x - 3.0) * (x - 3.0); };
        auto df = [](double x) { return 2.0 * (x - 3.0); };

        double xmin = 0.0;
        double fmin = dbrent(0.0, 3.0, 6.0, f, df, 1.0e-8, xmin);

        std::printf("Brent's method with derivatives (dbrent)\n");
        std::printf("=========================================\n\n");
        std::printf("Test 1: f(x) = (x-3)^2\n");
        std::printf("  x_min    = %16.10f  (expected 3.0)\n", xmin);
        std::printf("  f(x_min) = %16.10e  (expected 0.0)\n", fmin);
        std::printf("  error    = %16.10e\n\n", std::fabs(xmin - 3.0));

        // Test 2: f(x) = x^2 + sin(x), minimum near x = -0.4502
        auto f2  = [](double x) { return x * x + std::sin(x); };
        auto df2 = [](double x) { return 2.0 * x + std::cos(x); };
        fmin = dbrent(-2.0, -0.5, 1.0, f2, df2, 1.0e-10, xmin);

        std::printf("Test 2: f(x) = x^2 + sin(x)\n");
        std::printf("  x_min    = %16.10f\n", xmin);
        std::printf("  f(x_min) = %16.10f\n", fmin);
        std::printf("  f'(xmin) = %16.10e  (should be ~0)\n\n",
                    2.0 * xmin + std::cos(xmin));

        bool pass = std::fabs(xmin - 3.0) < 1.0e-7;
        dbrent(0.0, 3.0, 6.0, f, df, 1.0e-8, xmin);
        pass = std::fabs(xmin - 3.0) < 1.0e-7;
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
