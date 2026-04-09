// Driver for parallel MATAR CONVLV -- matches the Numerical Recipes convlv.dem
// test case. Convolves a rectangular pulse with a response kernel, then
// validates the FFT-based result against direct circular convolution.
//
// Usage: ./convlv [N]
//   N  -- signal length, must be power of 2 and >= 16 (default: 16)

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include <matar.h>
#include "convlv.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        int N = 16;
        if (argc > 1) N = std::atoi(argv[1]);
        if (N < 16 || (N & (N - 1)) != 0) {
            std::fprintf(stderr, "N must be a power of 2 and >= 16 (got %d)\n", N);
            MATAR_FINALIZE();
            return 1;
        }

        const int M = 9;

        // ---- host-side signal & response initialization ----
        DFMatrixKokkos<double> data_arr(N);
        DFMatrixKokkos<double> respns_orig(M);   // kept pristine for verification
        DFMatrixKokkos<double> resp(N);           // padded copy passed to convlv
        DFMatrixKokkos<double> ans(2 * N);

        for (int i = 1; i <= N; i++) {
            double v = 0.0;
            if (i >= N / 2 - N / 8 && i <= N / 2 + N / 8) v = 1.0;
            data_arr.host(i) = v;
        }

        for (int i = 1; i <= M; i++) {
            double v = 0.0;
            if (i > 2 && i < 7) v = 1.0;
            respns_orig.host(i) = v;
            resp.host(i)        = v;
        }
        for (int i = M + 1; i <= N; i++) {
            resp.host(i) = 0.0;
        }
        for (int i = 1; i <= 2 * N; i++) {
            ans.host(i) = 0.0;
        }

        data_arr.update_device();
        respns_orig.update_device();
        resp.update_device();
        ans.update_device();

        // ---- timed convolution ----
        auto t0 = std::chrono::high_resolution_clock::now();

        convlv(data_arr, N, resp, M, 1, ans);
        MATAR_FENCE();

        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

        ans.update_host();
        data_arr.update_host();
        respns_orig.update_host();

        // ---- verification: FFT result vs direct circular convolution ----
        std::printf("\n   I      CONVLV      Expected\n");
        double max_err = 0.0;
        for (int i = 1; i <= N; i++) {
            double cmp = 0.0;
            for (int j = 1; j <= M / 2; j++) {
                cmp += data_arr.host(((i - j - 1 + N) % N) + 1) * respns_orig.host(j + 1);
                cmp += data_arr.host(((i + j - 1)     % N) + 1) * respns_orig.host(M - j + 1);
            }
            cmp += data_arr.host(i) * respns_orig.host(1);

            double fft_val = ans.host(i);
            double err = std::fabs(fft_val - cmp);
            if (err > max_err) max_err = err;

            std::printf("  %2d   %12.6f %12.6f\n", i, fft_val, cmp);
        }

        std::printf("\n  Max |FFT - direct| = %.2e\n", max_err);
        std::printf("  CONVLV time: %.3f ms  (N = %d)\n\n", ms, N);
    }
    MATAR_FINALIZE();
    return 0;
}
