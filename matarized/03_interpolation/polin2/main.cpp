#include <cstdio>
#include <cmath>
#include <matar.h>
#include "polin2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        // Test with f(x1,x2) = x1^2 + x2^2
        constexpr int M = 4, N = 4;
        DFMatrixKokkos<double> x1a(M), x2a(N);
        DFMatrixKokkos<double> ya(M, N);

        for (int i = 1; i <= M; i++) x1a.host(i) = static_cast<double>(i);
        for (int j = 1; j <= N; j++) x2a.host(j) = static_cast<double>(j);

        for (int i = 1; i <= M; i++)
            for (int j = 1; j <= N; j++)
                ya.host(i, j) = x1a.host(i) * x1a.host(i) +
                                x2a.host(j) * x2a.host(j);

        double x1_test[] = {1.5, 2.5, 3.5};
        double x2_test[] = {1.5, 2.5, 3.5};
        constexpr int NTEST = 3;

        std::printf("2D polynomial interpolation of f(x1,x2) = x1^2 + x2^2\n\n");
        std::printf("%8s %8s %14s %14s %14s\n",
                    "x1", "x2", "polin2", "exact", "error");

        double max_err = 0.0;
        for (int t = 0; t < NTEST; t++) {
            double y, dy;
            polin2(x1a, x2a, ya, M, N, x1_test[t], x2_test[t], y, dy);
            double exact = x1_test[t] * x1_test[t] + x2_test[t] * x2_test[t];
            double err   = std::fabs(y - exact);
            if (err > max_err) max_err = err;
            std::printf("%8.2f %8.2f %14.8f %14.8f %14.2e\n",
                        x1_test[t], x2_test[t], y, exact, err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-8 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
