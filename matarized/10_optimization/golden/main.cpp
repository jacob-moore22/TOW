#include <cstdio>
#include <cmath>
#include <matar.h>
#include "golden.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Test 1: f(x) = (x - 3)^2, minimum at x = 3
        auto f1 = [](double x) { return (x - 3.0) * (x - 3.0); };
        double xmin = 0.0;
        double fmin = golden(0.0, 1.0, 6.0, f1, 1.0e-8, xmin);

        std::printf("Golden section search\n");
        std::printf("=====================\n\n");
        std::printf("Test 1: f(x) = (x-3)^2\n");
        std::printf("  x_min    = %16.10f  (expected 3.0)\n", xmin);
        std::printf("  f(x_min) = %16.10e  (expected 0.0)\n", fmin);
        std::printf("  error    = %16.10e\n\n", std::fabs(xmin - 3.0));

        // Test 2: f(x) = x^4 - 14x^3 + 60x^2 - 70x, minimum near x = 0.7808
        auto f2 = [](double x) {
            return x * x * x * x - 14.0 * x * x * x + 60.0 * x * x - 70.0 * x;
        };
        fmin = golden(0.0, 0.5, 2.0, f2, 1.0e-8, xmin);

        std::printf("Test 2: f(x) = x^4 - 14x^3 + 60x^2 - 70x\n");
        std::printf("  x_min    = %16.10f\n", xmin);
        std::printf("  f(x_min) = %16.10f\n\n", fmin);

        bool pass = std::fabs(xmin - 3.0) < 1.0e-6;
        // Rerun test 1 for pass/fail
        golden(0.0, 1.0, 6.0, f1, 1.0e-8, xmin);
        pass = std::fabs(xmin - 3.0) < 1.0e-6;
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
