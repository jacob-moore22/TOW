#include <cstdio>
#include <cmath>
#include <matar.h>
#include "pinvs.hpp"

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
                    c.host(i, j, k) = 0.0;

        s.host(1, 1) = 2.0;  s.host(1, 2) = 1.0;  s.host(1, 3) = 0.0;
        s.host(2, 1) = 1.0;  s.host(2, 2) = 3.0;  s.host(2, 3) = 1.0;
        s.host(3, 1) = 0.0;  s.host(3, 2) = 1.0;  s.host(3, 3) = 2.0;

        for (int i = 1; i <= NE; i++)
            s.host(i, NSJ) = static_cast<double>(i);

        std::printf("pinvs test: NE=%d\n", NE);
        pinvs(1, NE, 1, NSJ, 1, 1, c, NCI, NCJ, NCK, s, NSI, NSJ);

        std::printf("Reduced c(*,*,1) after pinvs:\n");
        for (int i = 1; i <= NE; i++) {
            for (int j = 1; j <= NCJ; j++)
                std::printf(" %8.4f", c.host(i, j, 1));
            std::printf("\n");
        }
        std::printf("Test PASSED (pinvs completed without error)\n");
    }
    MATAR_FINALIZE();
    return 0;
}
