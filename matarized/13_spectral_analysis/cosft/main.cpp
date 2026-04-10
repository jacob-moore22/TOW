// Driver for MATAR cosft -- discrete cosine transform.
// Tests forward transform of a pure cosine mode and round-trip recovery.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "cosft.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        int n = 16;
        if (argc > 1) n = std::atoi(argv[1]);
        if (n < 4 || (n & (n - 1)) != 0) {
            std::fprintf(stderr, "N must be a power of 2 and >= 4 (got %d)\n", n);
            MATAR_FINALIZE();
            return 1;
        }

        const double PI = 3.14159265358979;

        // --- Test 1: forward cosft of a known signal ---
        DFMatrixKokkos<double> y(n);
        for (int j = 1; j <= n; j++)
            y.host(j) = std::cos(PI * 2.0 * (j - 1) / n);
        y.update_device();

        cosft(y, n, 1);
        y.update_host();

        std::printf("Test 1: forward cosft of cos(2*pi*(j-1)/N), N=%d\n", n);
        for (int j = 1; j <= n; j++)
            std::printf("  %3d   %12.6f\n", j, y.host(j));

        // --- Test 2: round-trip (forward + inverse) ---
        DFMatrixKokkos<double> y2(n);
        DFMatrixKokkos<double> y2_orig(n);
        for (int j = 1; j <= n; j++) {
            double val = std::sin(2.0 * PI * j / n) + 0.3 * std::cos(6.0 * PI * j / n);
            y2.host(j) = val;
            y2_orig.host(j) = val;
        }
        y2.update_device();
        y2_orig.update_device();

        cosft(y2, n, 1);
        cosft(y2, n, -1);

        y2.update_host();
        y2_orig.update_host();

        double max_err = 0.0;
        double scale = n / 2.0;
        for (int j = 1; j <= n; j++) {
            double err = std::fabs(y2.host(j) / scale - y2_orig.host(j));
            if (err > max_err) max_err = err;
        }
        std::printf("\nTest 2: round-trip  max|cosft(-1, cosft(+1, y))/(N/2) - y| = %.2e\n\n",
                    max_err);
    }
    MATAR_FINALIZE();
    return 0;
}
