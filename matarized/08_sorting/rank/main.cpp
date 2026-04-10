#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "rank.hpp"
#include "../indexx/indexx.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 12;

        DFMatrixKokkos<double> arrin(N);
        DFMatrixKokkos<int>    indx(N);
        DFMatrixKokkos<int>    irank(N);

        std::srand(44444);
        std::printf("Original array:\n");
        for (int i = 1; i <= N; i++) {
            arrin.host(i) = static_cast<double>(std::rand() % 1000) / 10.0;
            std::printf("  arrin(%2d) = %8.2f\n", i, arrin.host(i));
        }

        indexx(N, arrin, indx);
        rank(N, indx, irank);

        std::printf("\nIndex and rank tables:\n");
        std::printf("%6s %10s %6s %6s\n", "j", "arrin(j)", "indx", "rank");
        for (int j = 1; j <= N; j++) {
            std::printf("%6d %10.2f %6d %6d\n",
                        j, arrin.host(j), indx.host(j), irank.host(j));
        }

        bool passed = true;
        for (int j = 1; j <= N; j++) {
            if (indx.host(irank.host(j)) != j) passed = false;
        }

        std::printf("\nRank-index consistency: %s\n", passed ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
