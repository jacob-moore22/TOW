#include <cstdio>
#include <cmath>
#include <matar.h>
#include "realft.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const int n = 16;
        DFMatrixKokkos<double> data(2 * n);

        for (int i = 1; i <= n; i++) {
            data.host(i) = std::sin(2.0 * M_PI * (i - 1) / n);
        }
        data.update_device();

        realft(data, n / 2, 1);
        MATAR_FENCE();
        data.update_host();

        std::printf("realft of sin(2*pi*t/N) with N=%d:\n", n);
        std::printf("  DC component (data[1]): %.6f (expect ~0)\n", data.host(1));
        std::printf("  Nyquist (data[2]):      %.6f (expect ~0)\n", data.host(2));

        double max_mag = 0.0;
        int peak_k = 0;
        for (int k = 1; k < n / 2; k++) {
            double re = data.host(2 * k + 1);
            double im = data.host(2 * k + 2);
            double mag = std::sqrt(re * re + im * im);
            if (mag > max_mag) { max_mag = mag; peak_k = k; }
        }
        std::printf("  Peak at k=%d, magnitude=%.6f (expect k=1, mag=%.1f)\n",
                    peak_k, max_mag, n / 2.0);

        realft(data, n / 2, -1);
        MATAR_FENCE();
        data.update_host();

        double max_err = 0.0;
        for (int i = 1; i <= n; i++) {
            double expected = std::sin(2.0 * M_PI * (i - 1) / n);
            double err = std::fabs(data.host(i) / (n / 2.0) - expected);
            if (err > max_err) max_err = err;
        }

        std::printf("\nRound-trip max error: %.2e\n", max_err);
        std::printf("Test %s\n", (peak_k == 1 && max_err < 1e-10) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
