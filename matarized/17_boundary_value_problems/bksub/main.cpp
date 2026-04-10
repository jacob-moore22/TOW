#include <cstdio>
#include <cmath>
#include <matar.h>
#include "bksub.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NE  = 3;
        constexpr int NB  = 1;
        constexpr int M   = 5;
        constexpr int NBF = NE - NB;
        constexpr int NCI = NE;
        constexpr int NCJ = NE - NB + 1;
        constexpr int NCK = M + 1;

        DFMatrixKokkos<double> c(NCI, NCJ, NCK);
        for (int i = 1; i <= NCI; i++)
            for (int j = 1; j <= NCJ; j++)
                for (int k = 1; k <= NCK; k++)
                    c.host(i, j, k) = 0.0;

        int jf = NCJ;
        for (int k = 1; k <= M; k++) {
            for (int i = 1; i <= NE; i++)
                c.host(i, jf, k) = static_cast<double>(i + k);
            for (int j = 1; j <= NBF; j++)
                for (int i = 1; i <= NE; i++)
                    c.host(i, j, k) = (i == j) ? 1.0 : 0.0;
        }
        for (int i = 1; i <= NE; i++)
            c.host(i, jf, M + 1) = static_cast<double>(i);

        std::printf("bksub test: NE=%d, NB=%d, M=%d\n", NE, NB, M);
        bksub(NE, NB, jf, 1, M, c, NCI, NCJ, NCK);

        std::printf("Solution column c(*,1,k) after bksub:\n");
        for (int k = 1; k <= M; k++) {
            std::printf("  k=%d:", k);
            for (int i = 1; i <= NE; i++)
                std::printf(" %8.4f", c.host(i, 1, k));
            std::printf("\n");
        }
        std::printf("Test PASSED (bksub completed without error)\n");
    }
    MATAR_FINALIZE();
    return 0;
}
