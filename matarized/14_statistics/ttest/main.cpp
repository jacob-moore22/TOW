#include <cstdio>
#include <cmath>
#include <matar.h>
#include "ttest.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N1 = 100, N2 = 100;
        DFMatrixKokkos<double> data1(N1), data2(N2);

        // Two samples from same distribution: data1 = 1..100, data2 = 1..100
        for (int i = 1; i <= N1; i++) data1.host(i) = static_cast<double>(i);
        for (int i = 1; i <= N2; i++) data2.host(i) = static_cast<double>(i);

        double t, prob;
        ttest(data1, N1, data2, N2, t, prob);

        std::printf("Student's t-test (identical samples):\n");
        std::printf("  t    = %14.6f  (expected ~ 0)\n", t);
        std::printf("  prob = %14.6f  (expected ~ 1)\n", prob);

        // Shifted sample
        for (int i = 1; i <= N2; i++) data2.host(i) = static_cast<double>(i) + 20.0;

        ttest(data1, N1, data2, N2, t, prob);
        std::printf("\nStudent's t-test (shifted by 20):\n");
        std::printf("  t    = %14.6f\n", t);
        std::printf("  prob = %14.10f  (expected very small)\n", prob);

        bool pass = fabs(t) > 2.0 && prob < 0.01;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
