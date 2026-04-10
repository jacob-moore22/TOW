#include <cstdio>
#include <cmath>
#include <matar.h>
#include "brent.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Test 1: f(x) = (x - 3)^2, minimum at x = 3
        auto f1 = [](double x) { return (x - 3.0) * (x - 3.0); };
        double xmin = 0.0;
        double fmin = brent(0.0, 3.0, 6.0, f1, 1.0e-8, xmin);

        std::printf("Brent's 1-D minimization\n");
        std::printf("========================\n\n");
        std::printf("Test 1: f(x) = (x-3)^2\n");
        std::printf("  x_min    = %16.10f  (expected 3.0)\n", xmin);
        std::printf("  f(x_min) = %16.10e  (expected 0.0)\n", fmin);
        std::printf("  error    = %16.10e\n\n", std::fabs(xmin - 3.0));

        // Test 2: f(x) = sin(x), minimum near x = -pi/2 in [-3, 0]
        auto f2 = [](double x) { return std::sin(x); };
        fmin = brent(-3.0, -1.5, 0.0, f2, 1.0e-8, xmin);

        double expected = -M_PI / 2.0;
        std::printf("Test 2: f(x) = sin(x), bracket [-3, 0]\n");
        std::printf("  x_min    = %16.10f  (expected %.10f)\n", xmin, expected);
        std::printf("  f(x_min) = %16.10f  (expected -1.0)\n", fmin);
        std::printf("  error    = %16.10e\n\n", std::fabs(xmin - expected));

        bool pass = std::fabs(xmin - expected) < 1.0e-7;
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
