#include <cstdio>
#include <cmath>
#include <matar.h>
#include "rtsafe.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("=== RTSAFE: Safe Newton-Raphson with bisection fallback ===\n\n");

        constexpr double TOL = 1.0e-14;

        auto fd1 = [](double x, double& f, double& df) {
            f  = x * x - 2.0;
            df = 2.0 * x;
        };
        double root1 = rtsafe(fd1, 1.0, 2.0, TOL);
        std::printf("x^2 - 2 = 0:   root = %.15f  (exact = %.15f)  err = %.2e\n",
                    root1, std::sqrt(2.0), std::fabs(root1 - std::sqrt(2.0)));

        auto fd2 = [](double x, double& f, double& df) {
            f  = std::sin(x);
            df = std::cos(x);
        };
        double root2 = rtsafe(fd2, 3.0, 4.0, TOL);
        std::printf("sin(x) = 0:    root = %.15f  (exact = %.15f)  err = %.2e\n",
                    root2, M_PI, std::fabs(root2 - M_PI));

        auto fd3 = [](double x, double& f, double& df) {
            f  = x * x * x - x - 1.0;
            df = 3.0 * x * x - 1.0;
        };
        double root3 = rtsafe(fd3, 1.0, 2.0, TOL);
        std::printf("x^3-x-1 = 0:   root = %.15f  f(root) = %.2e\n",
                    root3, std::fabs(root3 * root3 * root3 - root3 - 1.0));

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
