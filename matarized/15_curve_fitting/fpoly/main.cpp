#include <cstdio>
#include <cmath>
#include <matar.h>
#include "fpoly.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NP = 5;
        DFMatrixKokkos<double> p(NP);

        double x = 2.0;
        fpoly(x, p, NP);

        std::printf("Polynomial Basis Functions (FPOLY)\n");
        std::printf("==================================\n\n");
        std::printf("x = %.1f,  np = %d\n\n", x, NP);

        double expected[] = {1.0, 2.0, 4.0, 8.0, 16.0};
        double max_err = 0.0;

        for (int j = 1; j <= NP; j++) {
            double err = std::fabs(p.host(j) - expected[j - 1]);
            if (err > max_err) max_err = err;
            std::printf("  p(%d) = %12.6f  (expected = %6.1f, err = %.2e)\n",
                        j, p.host(j), expected[j - 1], err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-12 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
