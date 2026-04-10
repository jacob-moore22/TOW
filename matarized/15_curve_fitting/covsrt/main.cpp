#include <cstdio>
#include <cmath>
#include <matar.h>
#include "covsrt.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int MA   = 4;
        constexpr int MFIT = 3;

        DFMatrixKokkos<double> covar(MA, MA);
        DFMatrixKokkos<int> lista(MA);

        // Fit parameters 1, 3, 4 (parameter 2 is frozen)
        lista.host(1) = 1;
        lista.host(2) = 3;
        lista.host(3) = 4;
        lista.host(4) = 2;

        // Fill the top-left MFIT x MFIT block with a symmetric matrix
        double vals[3][3] = {{1.0, 0.5, 0.3},
                             {0.5, 2.0, 0.7},
                             {0.3, 0.7, 3.0}};
        for (int i = 1; i <= MFIT; i++)
            for (int j = 1; j <= MFIT; j++)
                covar.host(i, j) = vals[i - 1][j - 1];

        std::printf("Covariance Matrix Sort (COVSRT)\n");
        std::printf("===============================\n\n");

        std::printf("Before covsrt (MFIT x MFIT block):\n");
        for (int i = 1; i <= MFIT; i++) {
            std::printf("  ");
            for (int j = 1; j <= MFIT; j++)
                std::printf("%8.3f", covar.host(i, j));
            std::printf("\n");
        }
        std::printf("lista = {%d, %d, %d}, frozen = {%d}\n\n",
                    lista.host(1), lista.host(2), lista.host(3), lista.host(4));

        covsrt(covar, MA, lista, MFIT);

        std::printf("After covsrt (full MA x MA matrix):\n");
        for (int i = 1; i <= MA; i++) {
            std::printf("  ");
            for (int j = 1; j <= MA; j++)
                std::printf("%8.3f", covar.host(i, j));
            std::printf("\n");
        }

        // Verify: row/col 2 (frozen) should be all zeros
        double frozen_sum = 0.0;
        for (int j = 1; j <= MA; j++)
            frozen_sum += std::fabs(covar.host(2, j)) + std::fabs(covar.host(j, 2));

        // Verify symmetry
        double sym_err = 0.0;
        for (int i = 1; i <= MA; i++)
            for (int j = 1; j <= MA; j++)
                sym_err = std::fmax(sym_err, std::fabs(covar.host(i, j) - covar.host(j, i)));

        std::printf("\nFrozen row/col sum: %.6f (expected 0)\n", frozen_sum);
        std::printf("Symmetry error:     %.2e\n", sym_err);
        std::printf("Test %s\n", (frozen_sum < 1e-12 && sym_err < 1e-12) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
