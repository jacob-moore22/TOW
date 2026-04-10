#include <cstdio>
#include <cmath>
#include <matar.h>
#include "kendl1.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 20;
        DFMatrixKokkos<double> data1(N), data2(N);

        // Perfect concordance
        for (int i = 1; i <= N; i++) {
            data1.host(i) = static_cast<double>(i);
            data2.host(i) = static_cast<double>(i);
        }

        double tau, z, prob;
        kendl1(data1, data2, N, tau, z, prob);

        std::printf("Kendall's tau (perfect concordance):\n");
        std::printf("  tau  = %14.6f  (expected 1.0)\n", tau);
        std::printf("  z    = %14.6f\n", z);
        std::printf("  prob = %14.10f  (expected ~ 0)\n", prob);

        // Perfect discordance
        for (int i = 1; i <= N; i++)
            data2.host(i) = static_cast<double>(N + 1 - i);

        kendl1(data1, data2, N, tau, z, prob);

        std::printf("\nKendall's tau (perfect discordance):\n");
        std::printf("  tau  = %14.6f  (expected -1.0)\n", tau);
        std::printf("  z    = %14.6f\n", z);
        std::printf("  prob = %14.10f  (expected ~ 0)\n", prob);

        bool pass = fabs(tau + 1.0) < 1e-10;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
