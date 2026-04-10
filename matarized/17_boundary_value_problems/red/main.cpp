#include <cstdio>
#include <matar.h>
#include "red.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NE  = 3;
        constexpr int NCI = NE;
        constexpr int NCJ = 2 * NE + 1;
        constexpr int NCK = 10;
        constexpr int NSI = NE;
        constexpr int NSJ = 2 * NE + 1;

        DFMatrixKokkos<double> c(NCI, NCJ, NCK);
        DFMatrixKokkos<double> s(NSI, NSJ);

        for (int i = 1; i <= NCI; i++)
            for (int j = 1; j <= NCJ; j++)
                for (int k = 1; k <= NCK; k++)
                    c.host(i, j, k) = static_cast<double>((i + j + k) % 5);

        for (int i = 1; i <= NSI; i++)
            for (int j = 1; j <= NSJ; j++)
                s.host(i, j) = static_cast<double>(i * j);

        std::printf("red test: NE=%d\n", NE);
        std::printf("s before red:\n");
        for (int i = 1; i <= NSI; i++) {
            for (int j = 1; j <= NSJ; j++)
                std::printf(" %6.1f", s.host(i, j));
            std::printf("\n");
        }

        red(1, NE, 1, 2, 3, 4, NSJ, 1, 1, NCJ, 1, c, NCI, NCJ, NCK, s, NSI, NSJ);

        std::printf("s after red:\n");
        for (int i = 1; i <= NSI; i++) {
            for (int j = 1; j <= NSJ; j++)
                std::printf(" %6.1f", s.host(i, j));
            std::printf("\n");
        }
        std::printf("Test PASSED (red completed without error)\n");
    }
    MATAR_FINALIZE();
    return 0;
}
