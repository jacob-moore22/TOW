#include <cstdio>
#include <cmath>
#include <matar.h>
#include "twofft.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const int n = 16;
        DFMatrixKokkos<double> data1(n), data2(n);
        DFMatrixKokkos<double> fft1(2 * n), fft2(2 * n);

        for (int i = 1; i <= n; i++) {
            data1.host(i) = std::cos(2.0 * M_PI * (i - 1) / n);
            data2.host(i) = std::sin(2.0 * M_PI * 2.0 * (i - 1) / n);
        }
        data1.update_device();
        data2.update_device();

        twofft(data1, data2, fft1, fft2, n);
        MATAR_FENCE();
        fft1.update_host();
        fft2.update_host();

        std::printf("twofft of cos(2*pi*t/N) and sin(2*pi*2*t/N):\n");

        double max_mag1 = 0.0;
        int peak1 = 0;
        for (int k = 0; k < n; k++) {
            double re = fft1.host(2 * k + 1);
            double im = fft1.host(2 * k + 2);
            double mag = std::sqrt(re * re + im * im);
            if (mag > max_mag1) { max_mag1 = mag; peak1 = k; }
        }
        std::printf("  FFT1 peak at k=%d, magnitude=%.4f (expect k=1)\n", peak1, max_mag1);

        double max_mag2 = 0.0;
        int peak2 = 0;
        for (int k = 0; k < n; k++) {
            double re = fft2.host(2 * k + 1);
            double im = fft2.host(2 * k + 2);
            double mag = std::sqrt(re * re + im * im);
            if (mag > max_mag2) { max_mag2 = mag; peak2 = k; }
        }
        std::printf("  FFT2 peak at k=%d, magnitude=%.4f (expect k=2)\n", peak2, max_mag2);

        std::printf("\nTest %s\n", (peak1 == 1 && peak2 == 2) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
