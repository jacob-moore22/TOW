#include <cstdio>
#include <cmath>
#include <matar.h>
#include "crank.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 8;
        DFMatrixKokkos<double> w(N);

        // Sorted array with ties: {1, 2, 2, 4, 5, 5, 5, 8}
        double vals[] = {1.0, 2.0, 2.0, 4.0, 5.0, 5.0, 5.0, 8.0};
        for (int i = 1; i <= N; i++)
            w.host(i) = vals[i - 1];

        double s;
        crank(N, w, s);

        std::printf("Crank test (sorted with ties):\n");
        std::printf("  Ranks:");
        for (int i = 1; i <= N; i++)
            std::printf(" %.1f", w.host(i));
        std::printf("\n");
        std::printf("  S (tie correction) = %.1f\n", s);

        // Expected ranks: 1, 2.5, 2.5, 4, 6, 6, 6, 8
        // Tie group of 2: 2^3-2=6, group of 3: 3^3-3=24, total=30
        bool pass = w.host(1) == 1.0 && w.host(2) == 2.5 && w.host(3) == 2.5 &&
                    w.host(4) == 4.0 && w.host(5) == 6.0 && w.host(6) == 6.0 &&
                    w.host(7) == 6.0 && w.host(8) == 8.0 && fabs(s - 30.0) < 1e-10;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
