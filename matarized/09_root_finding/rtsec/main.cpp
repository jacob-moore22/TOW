#include <cstdio>
#include <cmath>
#include <matar.h>
#include "rtsec.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("=== RTSEC: Secant method ===\n\n");

        constexpr double TOL = 1.0e-12;

        auto f1 = [](double x) { return x * x - 2.0; };
        double root1 = rtsec(f1, 1.0, 2.0, TOL);
        std::printf("x^2 - 2 = 0:   root = %.15f  (exact = %.15f)  err = %.2e\n",
                    root1, std::sqrt(2.0), std::fabs(root1 - std::sqrt(2.0)));

        auto f2 = [](double x) { return std::sin(x); };
        double root2 = rtsec(f2, 3.0, 4.0, TOL);
        std::printf("sin(x) = 0:    root = %.15f  (exact = %.15f)  err = %.2e\n",
                    root2, M_PI, std::fabs(root2 - M_PI));

        auto f3 = [](double x) { return x * x * x - x - 1.0; };
        double root3 = rtsec(f3, 1.0, 2.0, TOL);
        std::printf("x^3-x-1 = 0:   root = %.15f  f(root) = %.2e\n",
                    root3, std::fabs(f3(root3)));

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
