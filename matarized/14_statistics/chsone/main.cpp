#include <cstdio>
#include <cmath>
#include <matar.h>
#include "chsone.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NBINS = 10;
        DFMatrixKokkos<double> bins(NBINS), ebins(NBINS);

        // Observed matches expected perfectly
        for (int i = 1; i <= NBINS; i++) {
            bins.host(i)  = 100.0;
            ebins.host(i) = 100.0;
        }

        double df, chsq, prob;
        chsone(bins, ebins, NBINS, 0, df, chsq, prob);

        std::printf("Chi-square one-sample (perfect match):\n");
        std::printf("  df   = %8.1f\n", df);
        std::printf("  chsq = %14.6f  (expected 0)\n", chsq);
        std::printf("  prob = %14.6f  (expected 1)\n", prob);

        // Perturbed bins
        bins.host(1) = 130.0;
        bins.host(2) = 70.0;

        chsone(bins, ebins, NBINS, 0, df, chsq, prob);

        std::printf("\nChi-square one-sample (perturbed):\n");
        std::printf("  df   = %8.1f\n", df);
        std::printf("  chsq = %14.6f\n", chsq);
        std::printf("  prob = %14.6f\n", prob);

        bool pass = chsq > 0.0 && prob < 1.0;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
