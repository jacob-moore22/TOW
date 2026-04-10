#include <cstdio>
#include <cmath>
#include <matar.h>
#include "pearsn.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 100;
        DFMatrixKokkos<double> x(N), y(N);

        // Perfect positive correlation: y = 2x + 3
        for (int i = 1; i <= N; i++) {
            x.host(i) = static_cast<double>(i);
            y.host(i) = 2.0 * i + 3.0;
        }

        double r, prob, z;
        pearsn(x, y, N, r, prob, z);

        std::printf("Pearson correlation (perfect linear y=2x+3):\n");
        std::printf("  r    = %14.10f  (expected 1.0)\n", r);
        std::printf("  prob = %14.10f  (expected ~ 0)\n", prob);
        std::printf("  z    = %14.6f\n", z);

        // Uncorrelated: x vs x^2 centered
        for (int i = 1; i <= N; i++) {
            double xi = i - 50.5;
            x.host(i) = xi;
            y.host(i) = xi * xi;
        }

        pearsn(x, y, N, r, prob, z);
        std::printf("\nPearson correlation (x vs x^2, symmetric => r~0):\n");
        std::printf("  r    = %14.10f  (expected ~ 0)\n", r);
        std::printf("  prob = %14.6f\n", prob);

        bool pass = fabs(r) < 0.01;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
