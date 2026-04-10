#include <cstdio>
#include <cmath>
#include <matar.h>
#include "tptest.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 100;
        DFMatrixKokkos<double> data1(N), data2(N);

        // Paired samples: data2 = data1 + small constant
        for (int i = 1; i <= N; i++) {
            data1.host(i) = static_cast<double>(i);
            data2.host(i) = static_cast<double>(i) + 0.01;
        }

        double t, prob;
        tptest(data1, data2, N, t, prob);

        std::printf("Paired t-test (shift = 0.01, should be barely significant):\n");
        std::printf("  t    = %14.6f\n", t);
        std::printf("  prob = %14.6f\n", prob);

        // Large shift
        for (int i = 1; i <= N; i++)
            data2.host(i) = static_cast<double>(i) + 100.0;

        tptest(data1, data2, N, t, prob);
        std::printf("\nPaired t-test (shift = 100):\n");
        std::printf("  t    = %14.6f\n", t);
        std::printf("  prob = %14.10f  (expected very small)\n", prob);

        bool pass = fabs(t) > 2.0 && prob < 0.01;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
