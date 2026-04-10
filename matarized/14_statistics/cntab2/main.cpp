#include <cstdio>
#include <cmath>
#include <matar.h>
#include "cntab2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NI = 2, NJ = 2;
        DFMatrixKokkos<double> nn(NI, NJ);

        // Independent table
        nn.host(1, 1) = 50.0; nn.host(1, 2) = 50.0;
        nn.host(2, 1) = 50.0; nn.host(2, 2) = 50.0;

        double h, hx, hy, hygx, hxgy, uygx, uxgy, uxy;
        cntab2(nn, NI, NJ, h, hx, hy, hygx, hxgy, uygx, uxgy, uxy);

        std::printf("Entropy measures (independent 2x2 table):\n");
        std::printf("  H    = %14.6f\n", h);
        std::printf("  HX   = %14.6f\n", hx);
        std::printf("  HY   = %14.6f\n", hy);
        std::printf("  HYGX = %14.6f\n", hygx);
        std::printf("  HXGY = %14.6f\n", hxgy);
        std::printf("  UYGX = %14.6f  (expected ~ 0)\n", uygx);
        std::printf("  UXGY = %14.6f  (expected ~ 0)\n", uxgy);
        std::printf("  UXY  = %14.6f  (expected ~ 0)\n", uxy);

        // Perfectly associated
        nn.host(1, 1) = 100.0; nn.host(1, 2) = 0.0;
        nn.host(2, 1) = 0.0;   nn.host(2, 2) = 100.0;

        cntab2(nn, NI, NJ, h, hx, hy, hygx, hxgy, uygx, uxgy, uxy);

        std::printf("\nEntropy measures (perfect association):\n");
        std::printf("  UYGX = %14.6f  (expected ~ 1)\n", uygx);
        std::printf("  UXGY = %14.6f  (expected ~ 1)\n", uxgy);
        std::printf("  UXY  = %14.6f  (expected ~ 1)\n", uxy);

        bool pass = fabs(uygx - 1.0) < 0.01 && fabs(uxy - 1.0) < 0.01;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
