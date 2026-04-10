#include <cstdio>
#include <cmath>
#include <matar.h>
#include "ddpoly.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // p(x) = x^3 - 2x + 1  ->  c = {1, -2, 0, 1}
        constexpr int NC = 4;
        constexpr int ND = 3;

        DFMatrixKokkos<double> c(NC);
        DFMatrixKokkos<double> pd(ND);

        c(1) = 1.0;  c(2) = -2.0;  c(3) = 0.0;  c(4) = 1.0;

        double x = 1.0;
        ddpoly(c, NC, x, pd, ND);

        std::printf("Polynomial p(x) = x^3 - 2x + 1 evaluated at x = %.1f:\n\n", x);
        std::printf("  p(1)   = %12.6f  (expected  0.0)\n", pd(1));
        std::printf("  p'(1)  = %12.6f  (expected  1.0)\n", pd(2));
        std::printf("  p''(1) = %12.6f  (expected  6.0)\n", pd(3));

        double err1 = std::fabs(pd(1) - 0.0);
        double err2 = std::fabs(pd(2) - 1.0);
        double err3 = std::fabs(pd(3) - 6.0);
        double max_err = std::fmax(err1, std::fmax(err2, err3));

        std::printf("\nMax absolute error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-12 ? "PASSED" : "FAILED");

        std::printf("\nMulti-point evaluation:\n");
        std::printf("%12s %18s %18s %18s\n", "X", "p(X)", "p'(X)", "p''(X)");
        for (int i = -5; i <= 5; i++) {
            x = 0.5 * i;
            ddpoly(c, NC, x, pd, ND);
            std::printf("%12.4f %18.10f %18.10f %18.10f\n", x, pd(1), pd(2), pd(3));
        }
    }
    MATAR_FINALIZE();
    return 0;
}
