// Driver for MATAR memcof -- maximum-entropy method LP coefficients.
// Generates a simple AR(2) process and recovers the coefficients.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <matar.h>
#include "memcof.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        const int N = 256;
        const int M = 4;

        // Generate a synthetic AR(2) signal: x(n) = 0.6*x(n-1) - 0.3*x(n-2) + noise
        DFMatrixKokkos<double> data(N);
        data.host(1) = 1.0;
        data.host(2) = 0.6;
        for (int j = 3; j <= N; j++) {
            double noise = 0.01 * std::sin(13.7 * j);
            data.host(j) = 0.6 * data.host(j - 1) - 0.3 * data.host(j - 2) + noise;
        }
        data.update_device();

        DFMatrixKokkos<double> cof(M), wk1(N), wk2(N), wkm(M);
        for (int j = 1; j <= M; j++) {
            cof.host(j) = 0.0;
            wkm.host(j) = 0.0;
        }
        for (int j = 1; j <= N; j++) {
            wk1.host(j) = 0.0;
            wk2.host(j) = 0.0;
        }
        cof.update_device();
        wk1.update_device();
        wk2.update_device();
        wkm.update_device();

        double pm = 0.0;
        memcof(data, N, M, pm, cof, wk1, wk2, wkm);

        cof.update_host();
        std::printf("MEM coefficients (AR(2) signal, M=%d poles):\n", M);
        std::printf("  Expected: a1 ~ 0.6, a2 ~ -0.3, others ~ 0\n");
        for (int j = 1; j <= M; j++)
            std::printf("  cof(%d) = %12.6f\n", j, cof.host(j));
        std::printf("  pm     = %12.6e\n\n", pm);
    }
    MATAR_FINALIZE();
    return 0;
}
