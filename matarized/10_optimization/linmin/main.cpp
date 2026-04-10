#include <cstdio>
#include <cmath>
#include <matar.h>
#include "linmin.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Line minimization (linmin)\n");
        std::printf("=========================\n\n");

        // Test 1: f(x,y) = (x-3)^2 + (y-5)^2, min at (3,5)
        auto sphere = [](const double* v) {
            return (v[0] - 3.0) * (v[0] - 3.0) + (v[1] - 5.0) * (v[1] - 5.0);
        };

        double p[2]  = {0.0, 0.0};
        double xi[2] = {3.0, 5.0};   // direction toward (3,5)
        double fret = 0.0;

        std::printf("Start:     p = (%.4f, %.4f)\n", p[0], p[1]);
        std::printf("Direction: xi = (%.4f, %.4f)\n\n", xi[0], xi[1]);

        linmin(p, xi, 2, fret, sphere);

        std::printf("After linmin:\n");
        std::printf("  p    = (%12.8f, %12.8f)  (expected 3, 5)\n", p[0], p[1]);
        std::printf("  fret = %16.10e  (expected 0)\n", fret);
        std::printf("  xi   = (%12.8f, %12.8f)\n\n", xi[0], xi[1]);

        double err = std::sqrt((p[0] - 3.0) * (p[0] - 3.0) +
                               (p[1] - 5.0) * (p[1] - 5.0));
        bool pass = err < 1.0e-4;
        std::printf("Distance from minimum: %.2e\n", err);
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
