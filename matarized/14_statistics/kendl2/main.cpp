#include <cstdio>
#include <cmath>
#include <matar.h>
#include "kendl2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NI = 3, NJ = 3;
        DFMatrixKokkos<double> tab(NI, NJ);

        // Strong positive association (diagonal dominance)
        tab.host(1, 1) = 20.0; tab.host(1, 2) = 1.0;  tab.host(1, 3) = 0.0;
        tab.host(2, 1) = 1.0;  tab.host(2, 2) = 20.0; tab.host(2, 3) = 1.0;
        tab.host(3, 1) = 0.0;  tab.host(3, 2) = 1.0;  tab.host(3, 3) = 20.0;

        double tau, z, prob;
        kendl2(tab, NI, NJ, tau, z, prob);

        std::printf("Kendall's tau (contingency table, diagonal dominance):\n");
        std::printf("  tau  = %14.6f  (expected positive)\n", tau);
        std::printf("  z    = %14.6f\n", z);
        std::printf("  prob = %14.10f  (expected small)\n", prob);

        // No association (uniform)
        for (int i = 1; i <= NI; i++)
            for (int j = 1; j <= NJ; j++)
                tab.host(i, j) = 10.0;

        kendl2(tab, NI, NJ, tau, z, prob);

        std::printf("\nKendall's tau (uniform table):\n");
        std::printf("  tau  = %14.6f  (expected ~ 0)\n", tau);
        std::printf("  z    = %14.6f\n", z);
        std::printf("  prob = %14.6f  (expected ~ 1)\n", prob);

        bool pass = fabs(tau) < 0.05;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
