#include <cstdio>
#include <cmath>
#include <matar.h>
#include "zbrac.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        std::printf("=== ZBRAC: Bracket a root ===\n\n");

        auto f1 = [](double x) { return x * x - 2.0; };
        auto f2 = [](double x) { return std::sin(x); };
        auto f3 = [](double x) { return x * x * x - x - 1.0; };

        double x1, x2;
        bool ok;

        x1 = 0.0; x2 = 1.0;
        ok = zbrac(f1, x1, x2);
        std::printf("f(x) = x^2 - 2:  bracket [%.6f, %.6f]  success=%d\n", x1, x2, ok);

        x1 = 2.0; x2 = 4.0;
        ok = zbrac(f2, x1, x2);
        std::printf("f(x) = sin(x):   bracket [%.6f, %.6f]  success=%d\n", x1, x2, ok);

        x1 = 1.0; x2 = 2.0;
        ok = zbrac(f3, x1, x2);
        std::printf("f(x) = x^3-x-1:  bracket [%.6f, %.6f]  success=%d\n", x1, x2, ok);

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
