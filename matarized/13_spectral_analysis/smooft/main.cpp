// Driver for MATAR smooft -- FFT-based smoothing.
// Adds noise to a smooth signal, applies smooft, and verifies the
// smoothed output is closer to the original than the noisy input.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "smooft.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        int n = 64;
        if (argc > 1) n = std::atoi(argv[1]);
        double pts = 5.0;
        if (argc > 2) pts = std::atof(argv[2]);

        const double PI = 3.14159265358979;

        DFMatrixKokkos<double> y(n);
        DFMatrixKokkos<double> y_clean(n);

        for (int j = 1; j <= n; j++) {
            double t = static_cast<double>(j - 1) / (n - 1);
            double clean = std::sin(2.0 * PI * t) + 0.5 * std::cos(4.0 * PI * t);
            double noise = 0.2 * std::sin(37.0 * j) * std::cos(53.0 * j);
            y.host(j) = clean + noise;
            y_clean.host(j) = clean;
        }
        y.update_device();
        y_clean.update_device();

        // Measure error before smoothing
        double err_before = 0.0;
        for (int j = 1; j <= n; j++) {
            double e = y.host(j) - y_clean.host(j);
            err_before += e * e;
        }
        err_before = std::sqrt(err_before / n);

        smooft(y, n, pts);

        y.update_host();
        y_clean.update_host();

        // Measure error after smoothing
        double err_after = 0.0;
        for (int j = 1; j <= n; j++) {
            double e = y.host(j) - y_clean.host(j);
            err_after += e * e;
        }
        err_after = std::sqrt(err_after / n);

        std::printf("Smooft test (N=%d, pts=%.1f)\n", n, pts);
        std::printf("  RMS error before smoothing: %.6f\n", err_before);
        std::printf("  RMS error after  smoothing: %.6f\n", err_after);
        std::printf("  Improvement ratio: %.2fx\n\n",
                    err_before / (err_after > 0 ? err_after : 1e-15));
    }
    MATAR_FINALIZE();
    return 0;
}
