#include <cstdio>
#include <cmath>
#include <matar.h>
#include "fgauss.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NA = 3;
        DFMatrixKokkos<double> a(NA);
        DFMatrixKokkos<double> dyda(NA);
        double y;

        // Single Gaussian: amplitude=5, center=10, width=3
        a.host(1) = 5.0;
        a.host(2) = 10.0;
        a.host(3) = 3.0;

        double x = 10.0;
        fgauss(x, a, y, dyda, NA);

        // At center x=10: arg=0, ex=1, y=5
        std::printf("Gaussian Basis Functions (FGAUSS)\n");
        std::printf("=================================\n\n");
        std::printf("Parameters: A=%.1f, center=%.1f, width=%.1f\n", a.host(1), a.host(2), a.host(3));
        std::printf("Evaluate at x = %.1f (center)\n\n", x);
        std::printf("  y      = %12.7f  (expected = 5.0)\n", y);
        std::printf("  dy/dA  = %12.7f  (expected = 1.0)\n", dyda.host(1));
        std::printf("  dy/dc  = %12.7f  (expected = 0.0)\n", dyda.host(2));
        std::printf("  dy/dw  = %12.7f  (expected = 0.0)\n", dyda.host(3));

        double err_y  = std::fabs(y - 5.0);
        double err_dA = std::fabs(dyda.host(1) - 1.0);
        double err_dc = std::fabs(dyda.host(2));
        double err_dw = std::fabs(dyda.host(3));

        // Also test at an off-center point
        x = 13.0;
        fgauss(x, a, y, dyda, NA);
        double arg = (13.0 - 10.0) / 3.0;
        double expected_y = 5.0 * std::exp(-arg * arg);
        double err_off = std::fabs(y - expected_y);

        std::printf("\nEvaluate at x = %.1f (off-center)\n", x);
        std::printf("  y      = %12.7f  (expected = %12.7f)\n", y, expected_y);
        std::printf("  dy/dA  = %12.7f\n", dyda.host(1));
        std::printf("  dy/dc  = %12.7f\n", dyda.host(2));
        std::printf("  dy/dw  = %12.7f\n", dyda.host(3));

        double max_err = std::fmax(std::fmax(err_y, err_dA),
                                   std::fmax(std::fmax(err_dc, err_dw), err_off));
        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
