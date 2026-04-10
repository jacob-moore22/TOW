#include <cstdio>
#include <cmath>
#include <matar.h>
#include "splie2.hpp"
#include "splin2.hpp"

using namespace mtr;

int main(int argc, char* argv[])
{
    MATAR_INITIALIZE(argc, argv);
    {
        constexpr int M = 8, N = 8;
        DFMatrixKokkos<double> x1a(M), x2a(N);
        DFMatrixKokkos<double> ya(M, N), y2a(M, N);

        for (int i = 1; i <= M; i++)
            x1a.host(i) = (i - 1) * M_PI / (M - 1);
        for (int j = 1; j <= N; j++)
            x2a.host(j) = (j - 1) * M_PI / (N - 1);

        for (int i = 1; i <= M; i++)
            for (int j = 1; j <= N; j++)
                ya.host(i, j) = std::sin(x1a.host(i)) * std::cos(x2a.host(j));

        splie2(x1a, x2a, ya, M, N, y2a);

        double x1_test[] = {0.5, 1.0, 1.5, 2.0, 2.5};
        double x2_test[] = {0.3, 0.8, 1.2, 1.8, 2.3};
        constexpr int NTEST = 5;

        std::printf("2D spline evaluation of f(x1,x2) = sin(x1)*cos(x2)\n\n");
        std::printf("%8s %8s %14s %14s %14s\n",
                    "x1", "x2", "splin2", "exact", "error");

        double max_err = 0.0;
        for (int t = 0; t < NTEST; t++) {
            double y;
            splin2(x1a, x2a, ya, y2a, M, N, x1_test[t], x2_test[t], y);
            double exact = std::sin(x1_test[t]) * std::cos(x2_test[t]);
            double err   = std::fabs(y - exact);
            if (err > max_err) max_err = err;
            std::printf("%8.4f %8.4f %14.8f %14.8f %14.2e\n",
                        x1_test[t], x2_test[t], y, exact, err);
        }

        std::printf("\nMax error: %.2e\n", max_err);
        std::printf("Test %s\n", max_err < 1e-2 ? "PASSED" : "FAILED");
    }
    MATAR_FINALIZE();
    return 0;
}
