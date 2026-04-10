#include <cstdio>
#include <cmath>
#include <matar.h>
#include "ftest.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N1 = 100, N2 = 100;
        DFMatrixKokkos<double> data1(N1), data2(N2);

        // Same variance: 1..100
        for (int i = 1; i <= N1; i++) data1.host(i) = static_cast<double>(i);
        for (int i = 1; i <= N2; i++) data2.host(i) = static_cast<double>(i);

        double f, prob;
        ftest(data1, N1, data2, N2, f, prob);

        std::printf("F-test (same variance):\n");
        std::printf("  F    = %14.6f  (expected ~ 1)\n", f);
        std::printf("  prob = %14.6f  (expected ~ 1)\n", prob);

        // Different variance: scale data2 by 5
        for (int i = 1; i <= N2; i++) data2.host(i) = static_cast<double>(i) * 5.0;

        ftest(data1, N1, data2, N2, f, prob);
        std::printf("\nF-test (variance ratio ~25):\n");
        std::printf("  F    = %14.6f  (expected ~ 25)\n", f);
        std::printf("  prob = %14.10f  (expected very small)\n", prob);

        bool pass = f > 20.0 && prob < 0.01;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
