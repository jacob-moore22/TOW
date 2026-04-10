#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "indexx.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 15;

        DFMatrixKokkos<double> arrin(N);
        DFMatrixKokkos<int>    indx(N);

        std::srand(33333);
        std::printf("Original array:\n");
        for (int i = 1; i <= N; i++) {
            arrin.host(i) = static_cast<double>(std::rand() % 1000) / 10.0;
            std::printf("  arrin(%2d) = %8.2f\n", i, arrin.host(i));
        }

        indexx(N, arrin, indx);

        std::printf("\nIndex array (indirect heapsort):\n");
        std::printf("%6s %6s %10s\n", "j", "indx", "arrin[indx]");
        bool passed = true;
        for (int j = 1; j <= N; j++) {
            std::printf("%6d %6d %10.2f\n", j, indx.host(j),
                        arrin.host(indx.host(j)));
            if (j > 1 && arrin.host(indx.host(j)) < arrin.host(indx.host(j - 1))) {
                passed = false;
            }
        }

        std::printf("\nTest %s\n", passed ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
