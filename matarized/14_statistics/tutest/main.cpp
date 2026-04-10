#include <cstdio>
#include <cmath>
#include <matar.h>
#include "tutest.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N1 = 100, N2 = 50;
        DFMatrixKokkos<double> data1(N1), data2(N2);

        for (int i = 1; i <= N1; i++) data1.host(i) = static_cast<double>(i);
        for (int i = 1; i <= N2; i++) data2.host(i) = static_cast<double>(i);

        double t, prob;
        tutest(data1, N1, data2, N2, t, prob);

        std::printf("Welch t-test (different sizes, overlapping data):\n");
        std::printf("  t    = %14.6f\n", t);
        std::printf("  prob = %14.6f\n", prob);

        // Shift data2 significantly
        for (int i = 1; i <= N2; i++) data2.host(i) = static_cast<double>(i) + 50.0;

        tutest(data1, N1, data2, N2, t, prob);
        std::printf("\nWelch t-test (shifted by 50):\n");
        std::printf("  t    = %14.6f\n", t);
        std::printf("  prob = %14.10f  (expected very small)\n", prob);

        bool pass = fabs(t) > 2.0 && prob < 0.01;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
