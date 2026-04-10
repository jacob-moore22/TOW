#include <cstdio>
#include <matar.h>
#include "hunt.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 7;
        DFMatrixKokkos<double> xx(N);

        double table[] = {0.0, 1.0, 3.0, 5.0, 8.0, 13.0, 21.0};
        for (int i = 1; i <= N; i++)
            xx.host(i) = table[i - 1];

        std::printf("Search with initial guess (HUNT) in ordered table\n");
        std::printf("Table: ");
        for (int i = 1; i <= N; i++) std::printf("%.0f ", xx.host(i));
        std::printf("\n\n");

        double x_test[]   = {4.0, 6.0, 12.0, 0.5};
        int    jlo_init[] = {2,   4,   3,    5};
        int    j_expect[] = {3,   4,   5,    1};
        constexpr int NTEST = 4;

        std::printf("%10s %10s %8s %10s %8s\n",
                    "x", "jlo_init", "jlo", "expected", "status");
        bool all_pass = true;
        for (int t = 0; t < NTEST; t++) {
            int jlo = jlo_init[t];
            hunt(xx, N, x_test[t], jlo);
            bool ok = (jlo == j_expect[t]);
            if (!ok) all_pass = false;
            std::printf("%10.1f %10d %8d %10d %8s\n",
                        x_test[t], jlo_init[t], jlo, j_expect[t],
                        ok ? "OK" : "FAIL");
        }
        std::printf("\nTest %s\n", all_pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
