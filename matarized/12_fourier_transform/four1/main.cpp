#include <cstdio>
#include <cmath>
#include <matar.h>
#include "four1.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const int nn = 8;
        DFMatrixKokkos<double> data(2 * nn);

        for (int i = 1; i <= 2 * nn; i++)
            data.host(i) = 0.0;
        data.host(1) = 1.0;
        data.update_device();

        four1(data, nn, 1);
        MATAR_FENCE();
        data.update_host();

        std::printf("FFT of unit impulse (should all be 1+0i):\n");
        double max_err = 0.0;
        for (int k = 1; k <= nn; k++) {
            double re = data.host(2 * k - 1);
            double im = data.host(2 * k);
            double err = std::fabs(re - 1.0) + std::fabs(im);
            if (err > max_err) max_err = err;
            std::printf("  k=%d: %10.6f + %10.6fi\n", k, re, im);
        }

        for (int i = 1; i <= 2 * nn; i++)
            data.host(i) = 0.0;
        data.host(1) = 1.0;
        data.update_device();

        four1(data, nn, 1);
        MATAR_FENCE();
        four1(data, nn, -1);
        MATAR_FENCE();
        data.update_host();

        double max_round = 0.0;
        for (int k = 1; k <= nn; k++) {
            double expect_re = (k == 1) ? (double)nn : 0.0;
            double expect_im = 0.0;
            double err = std::fabs(data.host(2 * k - 1) - expect_re) +
                         std::fabs(data.host(2 * k) - expect_im);
            if (err > max_round) max_round = err;
        }

        std::printf("\nMax impulse error: %.2e\n", max_err);
        std::printf("Max round-trip error: %.2e\n", max_round);
        std::printf("Test %s\n", (max_err < 1e-10 && max_round < 1e-10) ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
