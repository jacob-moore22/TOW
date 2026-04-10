#include <cstdio>
#include <cmath>
#include <matar.h>
#include "frprmn.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Conjugate gradient minimization (frprmn)\n");
        std::printf("=========================================\n\n");

        // Rosenbrock: f(x,y) = (1-x)^2 + 100(y-x^2)^2
        auto rosenbrock = [](const double* v) {
            double a = 1.0 - v[0];
            double b = v[1] - v[0] * v[0];
            return a * a + 100.0 * b * b;
        };
        auto grad_rosen = [](const double* v, double* df) {
            double x = v[0], y = v[1];
            df[0] = -2.0 * (1.0 - x) + 200.0 * (y - x * x) * (-2.0 * x);
            df[1] = 200.0 * (y - x * x);
        };

        double p[2] = {-1.0, 1.0};
        int iter = 0;
        double fret = 0.0;

        frprmn(p, 2, 1.0e-10, iter, fret, rosenbrock, grad_rosen);

        std::printf("Rosenbrock function: min at (1, 1)\n");
        std::printf("  p    = (%12.8f, %12.8f)\n", p[0], p[1]);
        std::printf("  fret = %16.10e\n", fret);
        std::printf("  iter = %d\n\n", iter);

        double err = std::sqrt((p[0] - 1.0) * (p[0] - 1.0) +
                               (p[1] - 1.0) * (p[1] - 1.0));

        // Test 2: f(x) = (x-3)^2, 1D
        auto quad = [](const double* v) { return (v[0] - 3.0) * (v[0] - 3.0); };
        auto dquad = [](const double* v, double* df) { df[0] = 2.0 * (v[0] - 3.0); };
        double p2[1] = {0.0};
        frprmn(p2, 1, 1.0e-10, iter, fret, quad, dquad);

        std::printf("f(x) = (x-3)^2: min at x=3\n");
        std::printf("  p    = %12.8f\n", p2[0]);
        std::printf("  fret = %16.10e\n\n", fret);

        bool pass = err < 1.0e-4 && std::fabs(p2[0] - 3.0) < 1.0e-6;
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
