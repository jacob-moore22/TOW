// Driver for MATAR fourn -- multi-dimensional FFT.
// Tests a 2D FFT of a delta function (should produce a flat spectrum)
// and round-trip recovery.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "fourn.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        int nx = 4, ny = 4;
        if (argc > 1) nx = std::atoi(argv[1]);
        if (argc > 2) ny = std::atoi(argv[2]);

        int ntot = nx * ny;
        int nn_dims[2] = {nx, ny};

        // --- Test 1: delta → flat spectrum ---
        DFMatrixKokkos<double> data(2 * ntot);
        for (int i = 1; i <= 2 * ntot; i++)
            data.host(i) = 0.0;
        data.host(1) = 1.0;  // delta at (0,0): real part
        data.update_device();

        fourn(data, nn_dims, 2, 1);

        data.update_host();
        std::printf("Test 1: 2D FFT of delta (%dx%d) -- all values should be 1+0i\n",
                    nx, ny);
        double max_err = 0.0;
        for (int i = 0; i < ntot; i++) {
            double re = data.host(2 * i + 1);
            double im = data.host(2 * i + 2);
            double err = std::fabs(re - 1.0) + std::fabs(im);
            if (err > max_err) max_err = err;
            if (ntot <= 16)
                std::printf("  (%d,%d): %8.4f + %8.4fi\n",
                            i / ny, i % ny, re, im);
        }
        std::printf("  Max |spectrum - 1| = %.2e\n", max_err);

        // --- Test 2: round-trip ---
        DFMatrixKokkos<double> data2(2 * ntot);
        DFMatrixKokkos<double> data2_orig(2 * ntot);
        for (int i = 1; i <= 2 * ntot; i++) {
            double val = std::sin(0.7 * i) + 0.3 * std::cos(1.3 * i);
            data2.host(i) = val;
            data2_orig.host(i) = val;
        }
        data2.update_device();
        data2_orig.update_device();

        fourn(data2, nn_dims, 2, 1);
        fourn(data2, nn_dims, 2, -1);

        data2.update_host();
        data2_orig.update_host();

        double max_rt = 0.0;
        double scale = ntot;
        for (int i = 1; i <= 2 * ntot; i++) {
            double err = std::fabs(data2.host(i) - scale * data2_orig.host(i));
            if (err > max_rt) max_rt = err;
        }
        std::printf("\nTest 2: round-trip  max|IFFT(FFT(x)) - N*x| = %.2e  (N=%d)\n\n",
                    max_rt, ntot);
    }
    MATAR_FINALIZE();
    return 0;
}
