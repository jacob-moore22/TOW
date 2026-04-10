#include <cstdio>
#include <cmath>
#include <matar.h>
#include "pzextr.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NV   = 2;
        constexpr int NUSE = 7;

        int nseq[] = {0, 2, 4, 6, 8, 12, 16, 24};

        double exact1 = M_PI;
        double exact2 = M_E;

        DFMatrixKokkos<double> yest(NV), yz(NV), dy(NV);

        std::printf("PZEXTR polynomial extrapolation test\n");
        std::printf("Synthetic data: exact = (pi, e), error = polynomial in h^2\n\n");
        std::printf("%6s %16s %16s %12s %12s\n",
                    "iest", "yz1", "yz2", "err1", "err2");

        for (int iest = 1; iest <= NUSE; iest++) {
            double h    = 1.0 / nseq[iest];
            double xest = h * h;

            yest.host(1) = exact1 + 0.5 * xest + 0.2 * xest * xest;
            yest.host(2) = exact2 + 0.3 * xest + 0.1 * xest * xest;

            pzextr(iest, xest, yest, yz, dy, NV, NUSE);

            double err1 = std::fabs(yz.host(1) - exact1);
            double err2 = std::fabs(yz.host(2) - exact2);
            std::printf("%6d %16.10f %16.10f %12.2e %12.2e\n",
                        iest, yz.host(1), yz.host(2), err1, err2);
        }

        double final_err = std::max(std::fabs(yz.host(1) - exact1),
                                    std::fabs(yz.host(2) - exact2));

        std::printf("\nFinal extrapolated values:\n");
        std::printf("  yz1 = %16.12f  (exact pi = %16.12f)\n", yz.host(1), exact1);
        std::printf("  yz2 = %16.12f  (exact e  = %16.12f)\n", yz.host(2), exact2);
        std::printf("  Max error: %.2e\n", final_err);
        std::printf("  Test %s\n", final_err < 1e-10 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
