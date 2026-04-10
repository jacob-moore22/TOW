// Driver for MATAR correl -- FFT-based cross-correlation.
// Correlates two signals and verifies against direct computation.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "correl.hpp"

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

        DFMatrixKokkos<double> data1(n), data2(n);
        DFMatrixKokkos<double> ans(2 * n);

        // Signal 1: rectangular pulse
        for (int j = 1; j <= n; j++) {
            data1.host(j) = (j >= n / 4 && j <= 3 * n / 4) ? 1.0 : 0.0;
            data2.host(j) = (j >= n / 4 + 2 && j <= 3 * n / 4 + 2) ? 1.0 : 0.0;
        }
        for (int j = 1; j <= 2 * n; j++)
            ans.host(j) = 0.0;

        data1.update_device();
        data2.update_device();
        ans.update_device();

        correl(data1, data2, n, ans);

        ans.update_host();
        data1.update_host();
        data2.update_host();

        // Direct circular cross-correlation for verification
        // correl computes R(i) = sum_j data1(j) * data2(j-i+1), i.e. negative-lag convention
        std::printf("  I   CORREL(FFT)    Direct\n");
        double max_err = 0.0;
        for (int i = 1; i <= n; i++) {
            double direct = 0.0;
            for (int j = 1; j <= n; j++) {
                int idx = ((j - 1 - (i - 1) + n) % n) + 1;
                direct += data1.host(j) * data2.host(idx);
            }
            double fft_val = ans.host(i);
            double err = std::fabs(fft_val - direct);
            if (err > max_err) max_err = err;
            std::printf(" %2d   %12.6f %12.6f\n", i, fft_val, direct);
        }
        std::printf("\n  Max |FFT - direct| = %.2e\n\n", max_err);
    }
    MATAR_FINALIZE();
    return 0;
}
