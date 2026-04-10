#include <cstdio>
#include <cmath>
#include <matar.h>
#include "ksone.hpp"

using namespace mtr;

// CDF for uniform distribution on [0,1]
double uniform_cdf(double x) { return x; }

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int N = 200;
        DFMatrixKokkos<double> data(N);

        // Uniform samples on [0,1]
        for (int i = 1; i <= N; i++)
            data.host(i) = (i - 0.5) / N;

        double d, prob;
        ksone(data, N, uniform_cdf, d, prob);

        std::printf("KS one-sample (uniform data vs uniform CDF):\n");
        std::printf("  D    = %14.6f  (expected small)\n", d);
        std::printf("  prob = %14.6f  (expected ~ 1)\n", prob);

        // Non-uniform: clustered at 0
        for (int i = 1; i <= N; i++)
            data.host(i) = ((i - 0.5) / N) * ((i - 0.5) / N);

        ksone(data, N, uniform_cdf, d, prob);

        std::printf("\nKS one-sample (quadratic data vs uniform CDF):\n");
        std::printf("  D    = %14.6f\n", d);
        std::printf("  prob = %14.10f  (expected very small)\n", prob);

        bool pass = d > 0.1 && prob < 0.01;
        std::printf("\nTest %s\n", pass ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
