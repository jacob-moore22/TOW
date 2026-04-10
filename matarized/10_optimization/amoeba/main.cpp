#include <cstdio>
#include <cmath>
#include <matar.h>
#include "amoeba.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("Nelder-Mead simplex (amoeba)\n");
        std::printf("============================\n\n");

        // Test 1: Rosenbrock f(x,y) = (1-x)^2 + 100(y-x^2)^2, min at (1,1)
        constexpr int ndim = 2;
        constexpr int mpts = ndim + 1;

        auto rosenbrock = [](const double* x) {
            double a = 1.0 - x[0];
            double b = x[1] - x[0] * x[0];
            return a * a + 100.0 * b * b;
        };

        // Initial simplex: three vertices in 2D
        double p[mpts * ndim] = {
            -1.0,  1.0,
             2.0,  0.0,
             0.0, -1.0
        };
        double y[mpts];
        for (int i = 0; i < mpts; i++)
            y[i] = rosenbrock(&p[i * ndim]);

        int iter = 0;
        amoeba(p, y, ndim, 1.0e-10, rosenbrock, iter);

        int ilo = 0;
        for (int i = 1; i < mpts; i++)
            if (y[i] < y[ilo]) ilo = i;

        std::printf("Test 1: Rosenbrock function\n");
        std::printf("  x_min    = (%12.8f, %12.8f)  (expected 1, 1)\n",
                    p[ilo * ndim], p[ilo * ndim + 1]);
        std::printf("  f(x_min) = %16.10e  (expected 0)\n", y[ilo]);
        std::printf("  iters    = %d\n\n", iter);

        double err = std::sqrt((p[ilo * ndim] - 1.0) * (p[ilo * ndim] - 1.0) +
                               (p[ilo * ndim + 1] - 1.0) * (p[ilo * ndim + 1] - 1.0));

        // Test 2: f(x) = (x-3)^2, 1D
        constexpr int ndim2 = 1;
        constexpr int mpts2 = 2;
        auto quad = [](const double* x) { return (x[0] - 3.0) * (x[0] - 3.0); };
        double p2[mpts2 * ndim2] = { 0.0, 6.0 };
        double y2[mpts2];
        for (int i = 0; i < mpts2; i++)
            y2[i] = quad(&p2[i * ndim2]);

        int iter2 = 0;
        amoeba(p2, y2, ndim2, 1.0e-10, quad, iter2);

        int ilo2 = (y2[0] < y2[1]) ? 0 : 1;
        std::printf("Test 2: f(x) = (x-3)^2\n");
        std::printf("  x_min    = %16.10f  (expected 3.0)\n", p2[ilo2]);
        std::printf("  f(x_min) = %16.10e\n", y2[ilo2]);
        std::printf("  iters    = %d\n\n", iter2);

        bool pass = err < 1.0e-4 && std::fabs(p2[ilo2] - 3.0) < 1.0e-4;
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
