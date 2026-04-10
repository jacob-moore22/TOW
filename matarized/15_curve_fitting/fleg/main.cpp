#include <cstdio>
#include <cmath>
#include <matar.h>
#include "fleg.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NL = 5;
        DFMatrixKokkos<double> pl(NL);

        double x = 0.5;
        fleg(x, pl, NL);

        // Known Legendre polynomial values at x=0.5:
        // P0=1, P1=0.5, P2=-0.125, P3=-0.4375, P4=-0.2890625
        double expected[] = {1.0, 0.5, -0.125, -0.4375, -0.2890625};

        std::printf("Legendre Polynomial Basis Functions (FLEG)\n");
        std::printf("==========================================\n\n");
        std::printf("x = %.1f,  nl = %d\n\n", x, NL);

        double max_err = 0.0;
        for (int j = 1; j <= NL; j++) {
            double err = std::fabs(pl.host(j) - expected[j - 1]);
            if (err > max_err) max_err = err;
            std::printf("  P_%d(%.1f) = %12.7f  (expected = %12.7f, err = %.2e)\n",
                        j - 1, x, pl.host(j), expected[j - 1], err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
