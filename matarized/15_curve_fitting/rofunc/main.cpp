#include <cstdio>
#include <cmath>
#include <matar.h>
#include "sort.hpp"
#include "rofunc.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NDATA = 20;
        DFMatrixKokkos<double> x(NDATA);
        DFMatrixKokkos<double> y(NDATA);
        DFMatrixKokkos<double> arr(NDATA);
        double aa, abdev;

        // Generate y = 2 + 3*x with small noise
        unsigned seed = 11111;
        for (int j = 1; j <= NDATA; j++) {
            x.host(j) = static_cast<double>(j);
            seed = seed * 1103515245 + 12345;
            double noise = 0.1 * (((seed >> 16) & 0x7FFF) / 16384.0 - 1.0);
            y.host(j) = 2.0 + 3.0 * x.host(j) + noise;
        }

        std::printf("Robust Line Fitting Objective (ROFUNC)\n");
        std::printf("======================================\n\n");
        std::printf("Data: y = 2 + 3*x + noise,  ndata = %d\n\n", NDATA);

        // Evaluate at true slope b=3 and a wrong slope b=1
        double f_true = rofunc(3.0, NDATA, x, y, arr, aa, abdev);
        double abdev_true = abdev;
        std::printf("At true slope b=3.0:\n");
        std::printf("  rofunc = %12.6f\n", f_true);
        std::printf("  aa     = %12.6f  (expected ~2)\n", aa);
        std::printf("  abdev  = %12.6f  (expected small)\n", abdev_true);

        double f_wrong = rofunc(1.0, NDATA, x, y, arr, aa, abdev);
        std::printf("\nAt wrong slope b=1.0:\n");
        std::printf("  rofunc = %12.6f  (expected large positive)\n", f_wrong);
        std::printf("  aa     = %12.6f\n", aa);
        std::printf("  abdev  = %12.6f\n", abdev);

        double f_high = rofunc(5.0, NDATA, x, y, arr, aa, abdev);
        std::printf("\nAt high slope b=5.0:\n");
        std::printf("  rofunc = %12.6f  (expected large negative)\n", f_high);
        std::printf("  aa     = %12.6f\n", aa);
        std::printf("  abdev  = %12.6f\n", abdev);

        // Sign change between b=1 and b=5 brackets the true slope
        bool sign_change = (f_wrong * f_high < 0.0);
        bool aa_near_2 = (std::fabs(aa - 2.0) < 2.0);
        bool abdev_small = (abdev_true < 1.0);
        std::printf("\nSign change between b=1 and b=5: %s\n", sign_change ? "yes" : "no");
        std::printf("Test %s\n", (sign_change && aa_near_2 && abdev_small) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
