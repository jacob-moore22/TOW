#include <cstdio>
#include <cmath>
#include <matar.h>
#include "chstwo.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NBINS = 10;
        DFMatrixKokkos<double> bins1(NBINS), bins2(NBINS);

        // Identical distributions
        for (int i = 1; i <= NBINS; i++) {
            bins1.host(i) = 50.0;
            bins2.host(i) = 50.0;
        }

        double df, chsq, prob;
        chstwo(bins1, bins2, NBINS, 0, df, chsq, prob);

        std::printf("Chi-square two-sample (identical):\n");
        std::printf("  df   = %8.1f\n", df);
        std::printf("  chsq = %14.6f  (expected 0)\n", chsq);
        std::printf("  prob = %14.6f  (expected 1)\n", prob);

        // Different distributions
        bins1.host(1) = 80.0;  bins2.host(1) = 20.0;
        bins1.host(2) = 20.0;  bins2.host(2) = 80.0;

        chstwo(bins1, bins2, NBINS, 0, df, chsq, prob);

        std::printf("\nChi-square two-sample (different):\n");
        std::printf("  df   = %8.1f\n", df);
        std::printf("  chsq = %14.6f\n", chsq);
        std::printf("  prob = %14.6f\n", prob);

        bool pass = chsq > 0.0 && prob < 1.0;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
