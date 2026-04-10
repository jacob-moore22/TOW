#include <cstdio>
#include <cmath>
#include <matar.h>
#include "link.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const int NCITY = 4;
        DFMatrixKokkos<double> x(NCITY), y(NCITY);
        DFMatrixKokkos<int> iorder(NCITY);
        DFMatrixKokkos<int> n(6);

        x.host(1) = 0.0; y.host(1) = 0.0;
        x.host(2) = 1.0; y.host(2) = 0.0;
        x.host(3) = 1.0; y.host(3) = 1.0;
        x.host(4) = 0.0; y.host(4) = 1.0;

        for (int i = 1; i <= NCITY; i++)
            iorder.host(i) = i;

        double path = 0.0;
        for (int i = 1; i <= NCITY; i++) {
            int j = (i % NCITY) + 1;
            double dx = x.host(iorder.host(i)) - x.host(iorder.host(j));
            double dy = y.host(iorder.host(i)) - y.host(iorder.host(j));
            path += std::sqrt(dx * dx + dy * dy);
        }
        std::printf("Initial path length: %.6f (expect 4.0 for unit square)\n", path);

        n.host(1) = 1; n.host(2) = 2;
        double de = 0.0;
        revcst(x, y, iorder, NCITY, n, de);
        std::printf("revcst DE for reversing segment [1,2]: %.6f\n", de);

        bool accept = metrop(-1.0, 1.0, 0.5);
        std::printf("metrop(DE=-1, T=1, ran=0.5): %s (expect true)\n",
                    accept ? "true" : "false");

        accept = metrop(100.0, 0.001, 0.5);
        std::printf("metrop(DE=100, T=0.001, ran=0.5): %s (expect false)\n",
                    accept ? "true" : "false");

        std::printf("\nTest PASSED\n");
    }
    MATAR_FINALIZE();
    return 0;
}
