#include <cstdio>
#include <matar.h>
#include "locate.hpp"

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

        double x_test[] = {-1.0, 0.5, 4.0, 8.0, 15.0, 25.0};
        int    j_expect[] = {0, 1, 3, 4, 6, 7};
        constexpr int NTEST = 6;

        std::printf("Binary search (LOCATE) in ordered table\n");
        std::printf("Table: ");
        for (int i = 1; i <= N; i++) std::printf("%.0f ", xx.host(i));
        std::printf("\n\n");

        std::printf("%10s %8s %10s %8s\n", "x", "j", "expected", "status");
        bool all_pass = true;
        for (int t = 0; t < NTEST; t++) {
            int j;
            locate(xx, N, x_test[t], j);
            bool ok = (j == j_expect[t]);
            if (!ok) all_pass = false;
            std::printf("%10.1f %8d %10d %8s\n",
                        x_test[t], j, j_expect[t], ok ? "OK" : "FAIL");
        }
        std::printf("\nTest %s\n", all_pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
