#include <cstdio>
#include <cmath>
#include <matar.h>
#include "sort.hpp"
#include "rofunc.hpp"
#include "medfit.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int NDATA = 30;
        DFMatrixKokkos<double> x(NDATA), y(NDATA);

        // Generate y = 2 + 3*x with noise, plus a few outliers
        unsigned seed = 11111;
        for (int i = 1; i <= NDATA; i++) {
            x.host(i) = static_cast<double>(i);
            seed = seed * 1103515245 + 12345;
            double noise = 0.3 * (((seed >> 16) & 0x7FFF) / 16384.0 - 1.0);
            y.host(i) = 2.0 + 3.0 * x.host(i) + noise;
        }

        // Add outliers
        y.host(5)  = 100.0;
        y.host(15) = -50.0;
        y.host(25) = 200.0;

        double a, b, abdev;
        medfit(x, y, NDATA, a, b, abdev);

        std::printf("Median-Based Robust Line Fitting (MEDFIT)\n");
        std::printf("=========================================\n\n");
        std::printf("Data: y = 2 + 3*x + noise, with 3 outliers\n\n");
        std::printf("Fitted parameters:\n");
        std::printf("  a (intercept) = %10.6f  (true = 2.0)\n", a);
        std::printf("  b (slope)     = %10.6f  (true = 3.0)\n", b);
        std::printf("  abdev         = %10.6f\n", abdev);

        double err_a = std::fabs(a - 2.0);
        double err_b = std::fabs(b - 3.0);

        std::printf("\nRecovery: |a-2| = %.4f, |b-3| = %.4f\n", err_a, err_b);
        bool pass = (err_a < 2.0) && (err_b < 0.5);
        std::printf("Test %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
