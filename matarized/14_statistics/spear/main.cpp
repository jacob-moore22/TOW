#include <cstdio>
#include <cmath>
#include <matar.h>
#include "spear.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 30;
        DFMatrixKokkos<double> data1(N), data2(N), wksp1(N), wksp2(N);

        // Perfect monotonic relationship
        for (int i = 1; i <= N; i++) {
            data1.host(i) = static_cast<double>(i);
            data2.host(i) = static_cast<double>(i) * 2.0 + 1.0;
        }

        double d, zd, probd, rs, probrs;
        spear(data1, data2, N, wksp1, wksp2, d, zd, probd, rs, probrs);

        std::printf("Spearman rank correlation (perfect monotonic):\n");
        std::printf("  D      = %14.6f  (expected 0)\n", d);
        std::printf("  zd     = %14.6f\n", zd);
        std::printf("  probd  = %14.6f\n", probd);
        std::printf("  rs     = %14.6f  (expected 1.0)\n", rs);
        std::printf("  probrs = %14.10f  (expected ~ 0)\n", probrs);

        // Reversed
        for (int i = 1; i <= N; i++) {
            data1.host(i) = static_cast<double>(i);
            data2.host(i) = static_cast<double>(N + 1 - i);
        }

        spear(data1, data2, N, wksp1, wksp2, d, zd, probd, rs, probrs);

        std::printf("\nSpearman rank correlation (reversed):\n");
        std::printf("  D      = %14.6f\n", d);
        std::printf("  rs     = %14.6f  (expected -1.0)\n", rs);
        std::printf("  probrs = %14.10f  (expected ~ 0)\n", probrs);

        bool pass = fabs(rs + 1.0) < 1e-6;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
