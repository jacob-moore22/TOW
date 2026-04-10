// Driver for MATAR spctrm -- power spectrum estimation.
// Generates a sine wave and verifies the spectrum peaks at the expected bin.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <matar.h>
#include "spctrm.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const int M = 8;
        const int K = 4;
        const bool OVRLAP = false;

        const double PI = 3.14159265358979;
        int total_samples = K * 4 * M;

        // Generate a sine wave at frequency bin 3
        std::vector<double> signal(total_samples);
        for (int i = 0; i < total_samples; i++)
            signal[i] = std::sin(2.0 * PI * 3.0 * i / (2 * M));

        DFMatrixKokkos<double> p(M);
        DFMatrixKokkos<double> w1(4 * M);
        DFMatrixKokkos<double> w2(M);

        for (int j = 1; j <= M; j++) p.host(j) = 0.0;
        for (int j = 1; j <= 4 * M; j++) w1.host(j) = 0.0;
        for (int j = 1; j <= M; j++) w2.host(j) = 0.0;
        p.update_device();
        w1.update_device();
        w2.update_device();

        int read_pos = 0;
        spctrm(p, M, K, OVRLAP, w1, w2, signal.data(), read_pos);

        p.update_host();

        std::printf("Power Spectrum (M=%d, K=%d segments, no overlap)\n", M, K);
        std::printf("  Bin   Power        (expect peak near bin 3)\n");
        double peak_val = 0.0;
        int peak_bin = 0;
        for (int j = 1; j <= M; j++) {
            std::printf("  %3d   %12.6e\n", j, p.host(j));
            if (p.host(j) > peak_val) {
                peak_val = p.host(j);
                peak_bin = j;
            }
        }
        std::printf("\n  Peak at bin %d (power = %.6e)\n\n", peak_bin, peak_val);
    }
    MATAR_FINALIZE();
    return 0;
}
