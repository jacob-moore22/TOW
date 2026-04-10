// Driver for MATAR sinft -- discrete sine transform.
// Verifies against a serial reference implementation and tests round-trip.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <matar.h>
#include "sinft.hpp"

using namespace mtr;

// Serial reference sinft (direct translation of Fortran)
static void sinft_ref(std::vector<double>& y, int n)
{
    double theta = 3.14159265358979 / static_cast<double>(n);
    double wr = 1.0, wi = 0.0;
    double wpr = -2.0 * std::sin(0.5 * theta) * std::sin(0.5 * theta);
    double wpi = std::sin(theta);
    y[1] = 0.0;
    int m = n / 2;
    for (int j = 1; j <= m; j++) {
        double wt = wr;
        wr = wr * wpr - wi * wpi + wr;
        wi = wi * wpr + wt * wpi + wi;
        double y1 = wi * (y[j + 1] + y[n - j + 1]);
        double y2 = 0.5 * (y[j + 1] - y[n - j + 1]);
        y[j + 1] = y1 + y2;
        y[n - j + 1] = y1 - y2;
    }
    // --- serial realft reference ---
    // Use MATAR realft via temporary device array
    DFMatrixKokkos<double> tmp(n);
    for (int j = 1; j <= n; j++) tmp.host(j) = y[j];
    tmp.update_device();
    realft(tmp, m, 1);
    tmp.update_host();
    for (int j = 1; j <= n; j++) y[j] = tmp.host(j);

    double sum = 0.0;
    y[1] = 0.5 * y[1];
    y[2] = 0.0;
    for (int j = 1; j <= n - 1; j += 2) {
        sum += y[j];
        y[j] = y[j + 1];
        y[j + 1] = sum;
    }
}

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
        int k_mode = 3;

        // --- Test 1: MATAR sinft vs serial reference ---
        DFMatrixKokkos<double> y(n);
        std::vector<double> y_ref(n + 1);
        for (int j = 1; j <= n; j++) {
            double val = std::sin(PI * k_mode * j / n);
            y.host(j) = val;
            y_ref[j] = val;
        }
        y.update_device();

        sinft(y, n);
        sinft_ref(y_ref, n);

        y.update_host();
        std::printf("sinft of sin(pi*%d*j/N), N=%d\n", k_mode, n);
        std::printf("   j    MATAR        Reference\n");
        double max_err = 0.0;
        for (int j = 1; j <= n; j++) {
            double err = std::fabs(y.host(j) - y_ref[j]);
            if (err > max_err) max_err = err;
            std::printf("  %2d  %12.6f  %12.6f\n", j, y.host(j), y_ref[j]);
        }
        std::printf("\n  Max |MATAR - ref| = %.2e\n", max_err);

        // --- Test 2: round-trip (y(1) must be 0 -- DST boundary condition) ---
        DFMatrixKokkos<double> y2(n);
        DFMatrixKokkos<double> y2_orig(n);
        y2.host(1) = 0.0;
        y2_orig.host(1) = 0.0;
        for (int j = 2; j <= n; j++) {
            double val = std::sin(2.0 * PI * (j - 1) / n);
            y2.host(j) = val;
            y2_orig.host(j) = val;
        }
        y2.update_device();
        y2_orig.update_device();

        sinft(y2, n);
        sinft(y2, n);

        y2.update_host();
        y2_orig.update_host();

        double max_rt = 0.0;
        double scale = n / 2.0;
        for (int j = 1; j <= n; j++) {
            double err = std::fabs(y2.host(j) - scale * y2_orig.host(j));
            if (err > max_rt) max_rt = err;
        }
        std::printf("  Round-trip  max|sinft(sinft(y)) - (N/2)*y| = %.2e\n\n", max_rt);
    }
    MATAR_FINALIZE();
    return 0;
}
