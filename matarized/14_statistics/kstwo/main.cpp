#include <cstdio>
#include <cmath>
#include <matar.h>
#include "kstwo.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N1 = 100, N2 = 100;
        DFMatrixKokkos<double> data1(N1), data2(N2);

        // Same distribution
        for (int i = 1; i <= N1; i++) data1.host(i) = (i - 0.5) / N1;
        for (int i = 1; i <= N2; i++) data2.host(i) = (i - 0.5) / N2;

        double d, prob;
        kstwo(data1, N1, data2, N2, d, prob);

        std::printf("KS two-sample (same distribution):\n");
        std::printf("  D    = %14.6f  (expected small)\n", d);
        std::printf("  prob = %14.6f  (expected ~ 1)\n", prob);

        // Different distributions
        for (int i = 1; i <= N1; i++) data1.host(i) = (i - 0.5) / N1;
        for (int i = 1; i <= N2; i++) data2.host(i) = ((i - 0.5) / N2) * ((i - 0.5) / N2);

        kstwo(data1, N1, data2, N2, d, prob);

        std::printf("\nKS two-sample (different distributions):\n");
        std::printf("  D    = %14.6f\n", d);
        std::printf("  prob = %14.10f  (expected very small)\n", prob);

        bool pass = d > 0.1 && prob < 0.01;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
